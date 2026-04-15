"""Tests for QPX core layer: DuckDBEngine, LazyQuery, QueryResult, Parquet I/O."""

import pandas as pd
import pyarrow as pa

from qpx.core.convert import QueryResult
from qpx.core.engine import DuckDBEngine, create_engine
from qpx.core.parquet_io import (
    parquet_row_count,
    read_parquet_metadata,
    read_parquet_schema,
)
from qpx.core.query import LazyQuery


def test_duckdb_engine(feature_parquet, tmp_path):
    """Engine: create, custom config, register parquet, context manager, params."""
    # Basic creation
    engine = DuckDBEngine()
    try:
        result = engine.execute("SELECT 42 AS answer")
        assert result.fetchone()[0] == 42
    finally:
        engine.close()

    # Custom memory and threads
    engine = DuckDBEngine(memory_limit="1GB", threads=2)
    try:
        result = engine.execute("SELECT current_setting('threads') AS t")
        assert int(result.fetchone()[0]) == 2
    finally:
        engine.close()

    # Register parquet creates view
    engine = DuckDBEngine()
    try:
        engine.register_parquet("feature", feature_parquet)
        count = engine.execute("SELECT COUNT(*) FROM feature").fetchone()[0]
        assert count == 3
    finally:
        engine.close()

    # Register replaces existing
    from qpx.writers import FeatureWriter
    from tests.conftest import make_feature_record

    path2 = tmp_path / "other.feature.parquet"
    with FeatureWriter(path2) as w:
        w.write_batch([make_feature_record()])

    engine = DuckDBEngine()
    try:
        engine.register_parquet("feature", feature_parquet)
        assert engine.execute("SELECT COUNT(*) FROM feature").fetchone()[0] == 3
        engine.register_parquet("feature", path2)
        assert engine.execute("SELECT COUNT(*) FROM feature").fetchone()[0] == 1
    finally:
        engine.close()

    # Context manager
    with DuckDBEngine() as engine:
        assert engine.execute("SELECT 1").fetchone()[0] == 1

    # Factory function
    engine = create_engine(memory_limit="2GB", threads=1)
    try:
        assert isinstance(engine, DuckDBEngine)
        assert engine.execute("SELECT 'ok'").fetchone()[0] == "ok"
    finally:
        engine.close()

    # Parameterized queries
    engine = DuckDBEngine()
    try:
        assert engine.execute("SELECT ? + ? AS sum", [10, 20]).fetchone()[0] == 30
    finally:
        engine.close()


def test_lazy_query(feature_parquet):
    """LazyQuery: count, filter, execute, distinct_values."""
    engine = DuckDBEngine()
    engine.register_parquet("feature", feature_parquet)
    try:
        q = LazyQuery(engine, "feature")
        assert q.count() == 3
        assert q.filter("is_decoy = false").count() == 2

        rows = q.execute().fetchall()
        assert len(rows) == 3

        runs = q.distinct_values("run_file_name")
        assert set(runs) == {"run_01", "run_02"}
    finally:
        engine.close()


def test_query_result(feature_parquet):
    """QueryResult: to_df, to_arrow, fetchone, fetchall, iter."""
    engine = DuckDBEngine()
    engine.register_parquet("feature", feature_parquet)
    try:
        # to_df
        raw = engine.execute("SELECT sequence, charge FROM feature")
        result = QueryResult(raw)
        df = result.to_df()
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 3
        assert "sequence" in df.columns

        # to_arrow
        raw = engine.execute("SELECT sequence, charge FROM feature")
        result = QueryResult(raw)
        table = result.to_arrow()
        assert isinstance(table, pa.Table)
        assert table.num_rows == 3

        # fetchone
        raw = engine.execute("SELECT COUNT(*) FROM feature")
        assert QueryResult(raw).fetchone()[0] == 3

        # fetchall
        raw = engine.execute("SELECT sequence FROM feature ORDER BY sequence")
        assert len(QueryResult(raw).fetchall()) == 3

        # iter
        raw = engine.execute("SELECT sequence FROM feature")
        assert len(list(QueryResult(raw))) == 3
    finally:
        engine.close()


def test_parquet_io(feature_parquet, dataset_parquet):
    """Parquet I/O: read_metadata, read_schema, row_count."""
    meta = read_parquet_metadata(feature_parquet)
    assert isinstance(meta, dict)
    assert meta["file_type"] == "feature_file"
    assert "qpx_version" in meta
    assert "creator" in meta
    assert "uuid" in meta
    assert len(meta["uuid"]) > 0

    schema = read_parquet_schema(feature_parquet)
    assert isinstance(schema, pa.Schema)
    assert "sequence" in schema.names

    assert parquet_row_count(feature_parquet) == 3
    assert parquet_row_count(dataset_parquet) == 1
