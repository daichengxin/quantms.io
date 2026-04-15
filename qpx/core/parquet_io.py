"""Parquet I/O utilities: metadata reading, schema inspection, row counts."""

from __future__ import annotations

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


def read_parquet_metadata(file_path: str | Path) -> dict[str, str]:
    """Read key-value metadata from a Parquet file footer."""
    pf = pq.ParquetFile(file_path)
    meta = pf.schema_arrow.metadata or {}
    return {k.decode(): v.decode() for k, v in meta.items()}


def read_parquet_schema(file_path: str | Path) -> pa.Schema:
    """Read just the Arrow schema from a Parquet file (no data)."""
    pf = pq.ParquetFile(file_path)
    return pf.schema_arrow


def parquet_row_count(file_path: str | Path) -> int:
    """Get row count without reading data."""
    pf = pq.ParquetFile(file_path)
    return pf.metadata.num_rows
