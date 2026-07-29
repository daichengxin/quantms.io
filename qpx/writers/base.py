"""BaseWriter — batched Parquet writer with schema validation."""

from __future__ import annotations

import datetime
import uuid
from pathlib import Path
from typing import TYPE_CHECKING

import pyarrow as pa
import pyarrow.parquet as pq

from qpx._version import __version__

if TYPE_CHECKING:
    import pandas as pd


# High-entropy float leaves that compress poorly under the default
# dictionary/PLAIN encoding (≈1.1x with ZSTD on continuous values). Parquet's
# BYTE_STREAM_SPLIT encoding reorders the bytes of IEEE floats so that ZSTD can
# exploit the per-byte regularity, measuring 15-28% smaller per column on real
# data. Paths use pyarrow's dotted leaf convention (``col``, ``a.list.element``,
# nested struct fields joined by ``.``). Only leaves present in a given schema
# are touched, so this list can name columns from any QPX view.
#
# NOTE: this is deliberately NOT applied to low-cardinality / dictionary-friendly
# float columns (e.g. ``additional_scores.score_value``,
# ``posterior_error_probability``), where BYTE_STREAM_SPLIT measured *larger*
# than the default dictionary encoding. Those keep their default encoding.
_BYTE_STREAM_SPLIT_LEAVES: tuple[str, ...] = (
    "predicted_rt",
    "rt",
    "rt_start",
    "rt_stop",
    "calculated_mz",
    "observed_mz",
    "intensities.list.element.intensity",
    "additional_intensities.list.element.intensities.list.element.intensity_value",
)

# ZSTD level 9 (vs the pyarrow default of 3) buys a few extra percent at modest
# write-cost; only codecs that accept a level get one.
_COMPRESSION_LEVELS: dict[str, int] = {"zstd": 9, "gzip": 9}


def _schema_leaf_paths(schema: pa.Schema) -> list[str]:
    """Return every leaf column path in pyarrow's dotted convention."""

    def _walk(field: pa.Field, prefix: str) -> list[str]:
        dtype = field.type
        if pa.types.is_struct(dtype):
            out: list[str] = []
            for i in range(dtype.num_fields):
                child = dtype.field(i)
                out.extend(_walk(child, f"{prefix}.{child.name}"))
            return out
        if pa.types.is_list(dtype) or pa.types.is_large_list(dtype):
            return _walk(dtype.value_field, f"{prefix}.list.element")
        return [prefix]

    leaves: list[str] = []
    for top in schema:
        leaves.extend(_walk(top, top.name))
    return leaves


class BaseWriter:
    """
    Batched Parquet writer with schema validation.

    Usage:
        writer = FeatureWriter("output.feature.parquet", creator="quantms")
        writer.write_batch(records)   # list[dict] or pd.DataFrame
        writer.write_table(table)     # pa.Table
        writer.close()
    """

    _schema_class = None  # Override in subclasses

    def __init__(
        self,
        path: str | Path,
        creator: str = "qpx",
        software_provider: str = "qpx",
        compression: str = "zstd",
        batch_size: int = 1_000_000,
        scan_format: str | None = None,
        extra_columns: list[str] | None = None,
    ):
        self._path = Path(path)
        self._compression = compression
        self._batch_size = batch_size
        self._buffer: list[dict] = []
        self._writer: pq.ParquetWriter | None = None
        self._extra_columns = extra_columns or []

        # Build Parquet footer metadata
        self._file_metadata = {
            b"qpx_version": __version__.encode(),
            b"creator": creator.encode(),
            b"file_type": self._schema_class._file_type.encode(),
            b"creation_date": datetime.datetime.now().isoformat().encode(),
            b"uuid": str(uuid.uuid4()).encode(),
            b"software_provider": software_provider.encode(),
            b"compression_format": compression.encode(),
        }
        if scan_format:
            self._file_metadata[b"scan_format"] = scan_format.encode()

    @property
    def arrow_schema(self) -> pa.Schema:
        if self._extra_columns:
            base = self._schema_class.get_extended_arrow_schema(self._extra_columns)
        else:
            base = self._schema_class.get_arrow_schema()
        return base.with_metadata(self._file_metadata)

    def write_batch(self, records: list[dict]):
        """Accumulate records and flush when batch_size is reached."""
        self._buffer.extend(records)
        while len(self._buffer) >= self._batch_size:
            batch = self._buffer[: self._batch_size]
            self._buffer = self._buffer[self._batch_size :]
            self._write_arrow_batch(batch)

    def write_table(self, table: pa.Table):
        """Write a complete Arrow table. Validates schema."""
        errors = self._schema_class.validate(table)
        if errors:
            raise ValueError("Schema validation failed:\n" + "\n".join(errors))
        self._ensure_writer()
        self._writer.write_table(table)

    def write_dataframe(self, df: "pd.DataFrame"):
        """Write a pandas DataFrame. Converts to Arrow and validates."""
        table = pa.Table.from_pandas(df, schema=self.arrow_schema, preserve_index=False)
        self.write_table(table)

    def close(self):
        """Flush remaining buffer and close the Parquet file."""
        if self._buffer:
            self._write_arrow_batch(self._buffer)
            self._buffer = []
        if self._writer:
            self._writer.close()
            self._writer = None

    def _write_arrow_batch(self, records: list[dict]):
        batch = pa.RecordBatch.from_pylist(records, schema=self.arrow_schema)
        # Validate non-nullable columns
        for i, field in enumerate(self.arrow_schema):
            if not field.nullable and batch.column(i).null_count > 0:
                raise ValueError(
                    f"Column '{field.name}' has {batch.column(i).null_count} null(s) but is marked as required in the schema"
                )
        self._ensure_writer()
        self._writer.write_batch(batch)

    def _ensure_writer(self):
        if self._writer is None:
            self._writer = pq.ParquetWriter(
                str(self._path),
                schema=self.arrow_schema,
                **self._parquet_write_options(),
            )

    def _parquet_write_options(self) -> dict:
        """Build ParquetWriter kwargs: codec level + selective BYTE_STREAM_SPLIT.

        BYTE_STREAM_SPLIT is applied only to the high-entropy float leaves
        present in this writer's schema; dictionary encoding is kept on every
        other column (pyarrow forbids dictionary + explicit encoding on the same
        column). All changes are encoding-only and lossless — the schema and
        values are untouched, so existing readers (pyarrow, DuckDB) read the
        output unchanged.
        """
        options: dict = {"compression": self._compression}
        level = _COMPRESSION_LEVELS.get(str(self._compression).lower())
        if level is not None:
            options["compression_level"] = level

        leaves = _schema_leaf_paths(self.arrow_schema)
        bss = {leaf: "BYTE_STREAM_SPLIT" for leaf in _BYTE_STREAM_SPLIT_LEAVES if leaf in leaves}
        if bss:
            options["column_encoding"] = bss
            # Dictionary stays on for everything except the BSS leaves.
            options["use_dictionary"] = [leaf for leaf in leaves if leaf not in bss]
            options["version"] = "2.6"  # BYTE_STREAM_SPLIT requires format v2.6+
        return options

    @staticmethod
    def write_partitioned(
        table: pa.Table,
        output_dir: str | Path,
        partition_cols: list[str] | None = None,
        compression: str = "zstd",
    ) -> Path:
        """Write Arrow table as Hive-partitioned Parquet.

        Parameters
        ----------
        table : pa.Table
            Data to write.
        output_dir : str | Path
            Root directory for partitioned output.
        partition_cols : list[str] | None
            Columns to partition by. Defaults to ["run_file_name"].
        compression : str
            Compression algorithm. Defaults to "zstd".

        Returns
        -------
        Path
            The output directory.
        """
        import pyarrow.dataset as ds

        output_dir = Path(output_dir)
        cols = partition_cols or ["run_file_name"]
        part_schema = pa.schema([table.schema.field(c) for c in cols])
        partitioning = ds.partitioning(part_schema, flavor="hive")
        file_options = ds.ParquetFileFormat().make_write_options(compression=compression)
        ds.write_dataset(
            table,
            str(output_dir),
            format="parquet",
            partitioning=partitioning,
            existing_data_behavior="overwrite_or_ignore",
            file_options=file_options,
        )
        return output_dir

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
