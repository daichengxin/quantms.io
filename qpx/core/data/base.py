"""BaseStructure — lazy DuckDB-backed data structure over a Parquet file."""

from __future__ import annotations

from pathlib import Path
from typing import Iterator

import pandas as pd
import pyarrow as pa

from qpx.core.convert import QueryResult
from qpx.core.data.schema import ValidationResult
from qpx.core.engine import DuckDBEngine, create_engine
from qpx.core.parquet_io import read_parquet_metadata
from qpx.core.query import LazyQuery, _escape_sql_string


class BaseStructure:
    """
    Lazy DuckDB-backed data structure over a Parquet file.

    Thin wrapper on core.query.LazyQuery — domain-specific methods only.
    All query building, materialization, and I/O lives in qpx.core.
    """

    _schema_class = None  # Override in subclasses (e.g., FeatureSchema)

    def __init__(self, engine: DuckDBEngine, table_name: str, file_path: Path):
        self._engine = engine
        self._table_name = table_name
        self._file_path = file_path
        self._query = LazyQuery(engine, table_name)

    @classmethod
    def from_file(cls, path: str | Path, **engine_kwargs) -> BaseStructure:
        """Open a standalone Parquet file (creates its own DuckDB engine)."""
        path = Path(path)
        engine = create_engine(**engine_kwargs)
        table_name = cls._schema_class._view_name
        engine.register_parquet(table_name, path)
        return cls(engine=engine, table_name=table_name, file_path=path)

    # --- Schema ---
    @classmethod
    def schema(cls) -> pa.Schema:
        return cls._schema_class.get_arrow_schema()

    # --- Lazy query building (delegates to core.query.LazyQuery) ---
    def select(self, *columns: str) -> BaseStructure:
        return self._with_query(self._query.select(*columns))

    def filter(self, condition: str) -> BaseStructure:
        return self._with_query(self._query.filter(condition))

    def limit(self, n: int) -> BaseStructure:
        return self._with_query(self._query.limit(n))

    def join(self, other: BaseStructure, on: str | None = None) -> BaseStructure:
        if on is None:
            on = self._auto_join_key(other)
        return self._with_query(self._query.join(other.query, on=on))

    def _join_list_membership(
        self,
        other: BaseStructure,
        list_column: str,
        value_column: str,
    ) -> BaseStructure:
        """Join by expanding one of this structure's list columns."""
        return self._with_query(
            self._query.join_list_membership(
                other.query,
                list_column=list_column,
                value_column=value_column,
            )
        )

    def order_by(self, *columns: str, desc: bool = False) -> BaseStructure:
        return self._with_query(self._query.order_by(*columns, desc=desc))

    # --- Materialization (delegates to core.convert.QueryResult) ---
    def to_df(self) -> pd.DataFrame:
        """Execute query and return pandas DataFrame."""
        return QueryResult(self._query.execute()).to_df()

    def to_arrow(self) -> pa.Table:
        """Execute query and return Arrow Table."""
        return QueryResult(self._query.execute()).to_arrow()

    def to_polars(self):
        """Execute query and return Polars DataFrame."""
        return QueryResult(self._query.execute()).to_polars()

    def count(self) -> int:
        return self._query.count()

    @property
    def query(self) -> LazyQuery:
        """Return the immutable query backing this structure."""
        return self._query

    # --- Batch iteration (for sample-level processing) ---
    def iter_batches(self, partition_by: str, batch_size: int = 20) -> Iterator[tuple[list[str], pd.DataFrame]]:
        """
        Iterate over data in batches, partitioned by a column.

        Yields (partition_values, DataFrame) tuples. Useful for
        memory-bounded processing of large datasets.
        """
        all_values = self._query.distinct_values(partition_by)
        for i in range(0, len(all_values), batch_size):
            batch_values = all_values[i : i + batch_size]
            placeholders = ", ".join(f"'{_escape_sql_string(str(v))}'" for v in batch_values)
            filtered = self._query.filter(f"{partition_by} IN ({placeholders})")
            df = QueryResult(filtered.execute()).to_df()
            yield batch_values, df

    # --- Validation ---
    def validate(self, *, strict: bool = False) -> ValidationResult:
        """Validate this structure's data against its schema.

        ``strict=True`` (used by ``qpxc validate``) promotes duplicate-PK and
        required-null issues to errors; the default is lenient.
        """
        table = self.to_arrow()
        return self._schema_class.validate_full(table, strict=strict)

    # --- Parquet metadata (delegates to core.parquet_io) ---
    @property
    def file_metadata(self) -> dict[str, str]:
        return read_parquet_metadata(self._file_path)

    # --- Internals ---
    def _with_query(self, new_query: LazyQuery) -> BaseStructure:
        """Return a clone with a different LazyQuery (immutable pattern)."""
        clone = self.__class__.__new__(self.__class__)
        clone._engine = self._engine
        clone._table_name = self._table_name
        clone._file_path = self._file_path
        clone._query = new_query
        return clone

    def _auto_join_key(self, other: BaseStructure) -> str:
        my_cols = set(self.schema().names)
        their_cols = set(other.schema().names)
        common = (my_cols & their_cols) - {"is_decoy"}
        if "run_file_name" in common:
            return "run_file_name"
        if "sample_accession" in common:
            return "sample_accession"
        if common:
            return next(iter(common))
        raise ValueError("No common column found for auto-join. Specify 'on' explicitly.")

    def __repr__(self):
        return f"{self.__class__.__name__}('{self._file_path}', rows={self.count()})"

    def __len__(self):
        return self.count()
