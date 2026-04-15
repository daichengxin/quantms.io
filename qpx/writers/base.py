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
                compression=self._compression,
            )

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
