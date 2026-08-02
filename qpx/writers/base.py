"""BaseWriter — batched Parquet writer with schema validation."""

from __future__ import annotations

import datetime
import uuid
from pathlib import Path
from typing import TYPE_CHECKING

import pyarrow as pa
import pyarrow.parquet as pq

from qpx._version import __version__
from qpx.version import QPX_SPEC_VERSION

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
#
# PEP was a candidate for float32 + BYTE_STREAM_SPLIT (audit projected ~52%),
# but empirical validation over the real qpx corpus found genuine PEP values as
# small as ~1e-300..1e-318 (DIA/quantms feature+psm files): these underflow to
# zero in float32 and would corrupt any downstream -log10(PEP)/q-value use. The
# float32 narrowing was therefore held (PEP stays float64), so it is not split.
#
# The intensity leaves (``intensities.list.element.intensity`` and the nested
# ``additional_intensities…intensity_value``) are intentionally NOT split.
# BYTE_STREAM_SPLIT helps continuous LFQ intensities (~-20%) but slightly *hurts*
# TMT reporter-ion columns (+4-7%), which dominate TMT feature files (~64% of
# the file). Rather than probe the data at write time, we simply omit them here:
# they fall back to dictionary encoding, which is safe everywhere (never worse
# than PLAIN) and avoids the TMT regression. RT/mz leaves (a proven win) stay.
_BYTE_STREAM_SPLIT_LEAVES: tuple[str, ...] = (
    "predicted_rt",
    "rt",
    "rt_start",
    "rt_stop",
    "calculated_mz",
    "observed_mz",
)

# ZSTD level 9 (vs the pyarrow default of 3) buys a few extra percent at modest
# write-cost; only codecs that accept a level get one.
_COMPRESSION_LEVELS: dict[str, int] = {"zstd": 9, "gzip": 9}

# Row-group sizing for partitioned (dataset) writes. Left unset, pyarrow's
# dataset writer flushes a group per input batch, which for streamed msnet
# conversion fragments into ~80-row groups and collapses ZSTD (1.12x vs 1.75x).
# A 200k-1M window matches the single-file writer's effective grouping.
_PARTITIONED_MIN_ROWS_PER_GROUP: int = 200_000
_PARTITIONED_MAX_ROWS_PER_GROUP: int = 1_000_000


def _stamp_footer_metadata(
    table: pa.Table,
    compression: str,
    file_metadata: dict | None = None,
) -> pa.Table:
    """Return *table* with qpx footer KV metadata on its schema.

    Layering: any metadata already on the table schema, then defaults that a
    caller may override (``creation_date``/``uuid``), then *file_metadata* (e.g.
    a writer's ``_file_metadata`` carrying ``file_type``/``uuid``), and finally
    the **canonical identity** (``qpx_version``/``writer_version``/
    ``compression_format``) which always wins so a stale caller value cannot
    mis-stamp the file. Used so partitioned datasets carry the same footer
    identity as single files.
    """
    merged: dict[bytes, bytes] = dict(table.schema.metadata or {})
    # Overridable defaults.
    merged.update(
        {
            b"creation_date": datetime.datetime.now().isoformat().encode(),
            b"uuid": str(uuid.uuid4()).encode(),
        }
    )
    if file_metadata:
        merged.update(
            {
                (k if isinstance(k, bytes) else str(k).encode()): (v if isinstance(v, bytes) else str(v).encode())
                for k, v in file_metadata.items()
            }
        )
    # Canonical identity — reserved, must not be overridden by caller metadata.
    merged.update(
        {
            b"qpx_version": QPX_SPEC_VERSION.encode(),
            b"writer_version": __version__.encode(),
            b"compression_format": str(compression).encode(),
        }
    )
    return table.replace_schema_metadata(merged)


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


def parquet_write_options(schema: pa.Schema, compression: str) -> dict:
    """Build QPX ParquetWriter options for *schema* and *compression*."""
    options: dict = {"compression": compression}
    level = _COMPRESSION_LEVELS.get(str(compression).lower())
    if level is not None:
        options["compression_level"] = level

    leaves = _schema_leaf_paths(schema)
    bss = {leaf: "BYTE_STREAM_SPLIT" for leaf in _BYTE_STREAM_SPLIT_LEAVES if leaf in leaves}
    # Always stamp the explicit v2.6 format and dictionary-encode non-BSS leaves so
    # every artifact (feature/psm and pg alike) uses the same encoding, rather than
    # falling back to pyarrow defaults when no byte-stream-split leaf is present.
    options["version"] = "2.6"
    options["use_dictionary"] = [leaf for leaf in leaves if leaf not in bss]
    if bss:
        options["column_encoding"] = bss
    return options


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

        # Build Parquet footer metadata. Two distinct versions are stamped:
        #   qpx_version    — the on-disk *specification* version (QPX_SPEC_VERSION),
        #                    what readers check for format compatibility.
        #   writer_version — the *package* build (__version__, hatch-vcs), for
        #                    provenance/debugging only.
        self._file_metadata = {
            b"qpx_version": QPX_SPEC_VERSION.encode(),
            b"writer_version": __version__.encode(),
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
        return parquet_write_options(self.arrow_schema, self._compression)

    @staticmethod
    def write_partitioned(
        table: pa.Table,
        output_dir: str | Path,
        partition_cols: list[str] | None = None,
        compression: str = "zstd",
        file_metadata: dict | None = None,
    ) -> Path:
        """Write Arrow table as Hive-partitioned Parquet.

        Encoding parity with the single-file writers: the same
        :func:`parquet_write_options` path is used, so partitioned output gets
        the tuned ZSTD level, selective BYTE_STREAM_SPLIT and dictionary
        encoding, Parquet v2.6, plus a sane row-group size (so compression is
        not defeated by ~80-row groups). The qpx footer KV metadata is stamped
        onto every file so partitioned datasets keep ``qpx_version`` /
        ``file_type`` / ``uuid`` like single files do.

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
        file_metadata : dict | None
            Extra footer KV metadata (bytes->bytes) to stamp, e.g. a writer's
            ``_file_metadata`` carrying ``file_type`` / ``uuid``. Merged on top
            of any metadata already on the table schema and a minimal qpx footer.

        Returns
        -------
        Path
            The output directory.
        """
        import pyarrow.dataset as ds

        output_dir = Path(output_dir)
        cols = partition_cols or ["run_file_name"]

        # Validate the partition columns up front so the default path (or an
        # unsupported key) fails with an actionable message instead of a raw
        # KeyError deep in pyarrow. Since QPX 1.1 the pg view is keyed on
        # grouped_runs (list<string>) with no scalar run_file_name, and Hive
        # partitioning encodes each value as a directory name — it cannot
        # represent a list — so the pg view is not Hive-partitionable.
        available = set(table.schema.names)
        missing = [c for c in cols if c not in available]
        if missing:
            if "grouped_runs" in available and "run_file_name" in missing:
                raise ValueError(
                    "The protein-group (pg) view is not Hive-partitionable: it is keyed on "
                    "'grouped_runs' (list<string>, since QPX 1.1) and has no scalar "
                    "'run_file_name' column. Hive partitioning encodes each value as a "
                    "directory name and cannot represent a list. Write pg as a single "
                    "Parquet file (write_table), or partition the psm/feature views, which "
                    "keep 'run_file_name'."
                )
            raise ValueError(f"Cannot partition on missing column(s) {missing}; available columns: {sorted(available)}")
        for col in cols:
            col_type = table.schema.field(col).type
            if pa.types.is_list(col_type) or pa.types.is_large_list(col_type) or pa.types.is_struct(col_type):
                raise ValueError(
                    f"Cannot Hive-partition on non-scalar column '{col}' (type {col_type}); partition "
                    "columns must be scalar. The pg view's 'grouped_runs' is a list, so pg is not "
                    "partitionable — write it as a single file."
                )

        table = _stamp_footer_metadata(table, compression, file_metadata)

        part_schema = pa.schema([table.schema.field(c) for c in cols])
        partitioning = ds.partitioning(part_schema, flavor="hive")

        opts = parquet_write_options(table.schema, compression)
        file_options = ds.ParquetFileFormat().make_write_options(**opts)
        ds.write_dataset(
            table,
            str(output_dir),
            format="parquet",
            partitioning=partitioning,
            existing_data_behavior="overwrite_or_ignore",
            file_options=file_options,
            min_rows_per_group=_PARTITIONED_MIN_ROWS_PER_GROUP,
            max_rows_per_group=_PARTITIONED_MAX_ROWS_PER_GROUP,
        )
        return output_dir

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
