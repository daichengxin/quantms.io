"""Tests for enhanced views: BaseView caching, new views, Dataset view properties."""

from qpx.dataset import Dataset
from qpx.views.base import BaseView, CachedQueryResult


class TestBaseViewCaching:
    """Test BaseView caching and materialization."""

    def test_execute_cached_returns_cached_query_result(self, dataset_dir):
        with Dataset(dataset_dir) as ds:
            view = BaseView(ds)
            result = view._execute_cached("test_key", "SELECT COUNT(*) AS cnt FROM feature")
            assert isinstance(result, CachedQueryResult)

    def test_cache_returns_same_object(self, dataset_dir):
        with Dataset(dataset_dir) as ds:
            view = BaseView(ds)
            sql = "SELECT COUNT(*) AS cnt FROM feature"
            r1 = view._execute_cached("test_key", sql)
            r2 = view._execute_cached("test_key", sql)
            assert r1 is r2

    def test_invalidate_cache(self, dataset_dir):
        with Dataset(dataset_dir) as ds:
            view = BaseView(ds)
            sql = "SELECT COUNT(*) AS cnt FROM feature"
            r1 = view._execute_cached("test_key", sql)
            view.invalidate_cache()
            r2 = view._execute_cached("test_key", sql)
            assert r1 is not r2

    def test_to_parquet(self, dataset_dir, tmp_path):
        with Dataset(dataset_dir) as ds:
            view = BaseView(ds)
            result = view._execute_cached("cnt", "SELECT COUNT(*) AS cnt FROM feature")
            out = tmp_path / "test_out.parquet"
            view.to_parquet(result, out)
            assert out.exists()

    def test_cached_result_to_df(self, dataset_dir):
        with Dataset(dataset_dir) as ds:
            view = BaseView(ds)
            result = view._execute_cached("cnt", "SELECT COUNT(*) AS cnt FROM feature")
            df = result.to_df()
            assert len(df) == 1
            assert df.iloc[0]["cnt"] == 3

    def test_cached_result_to_arrow(self, dataset_dir):
        import pyarrow as pa

        with Dataset(dataset_dir) as ds:
            view = BaseView(ds)
            result = view._execute_cached("cnt", "SELECT COUNT(*) AS cnt FROM feature")
            table = result.to_arrow()
            assert isinstance(table, pa.Table)
            assert table.num_rows == 1


class TestRunSummaryView:
    """Test RunSummaryView."""

    def test_summary_returns_cached_result(self, dataset_dir):
        from qpx.views.api import RunSummaryView

        with Dataset(dataset_dir) as ds:
            view = RunSummaryView(ds)
            result = view.summary()
            assert isinstance(result, CachedQueryResult)

    def test_summary_has_expected_columns(self, dataset_dir):
        from qpx.views.api import RunSummaryView

        with Dataset(dataset_dir) as ds:
            view = RunSummaryView(ds)
            df = view.summary().to_df()
            assert "run_file_name" in df.columns
            assert "n_peptides" in df.columns
            assert "n_proteins" in df.columns
            assert "n_features" in df.columns

    def test_summary_aggregates_per_run(self, dataset_dir):
        from qpx.views.api import RunSummaryView

        with Dataset(dataset_dir) as ds:
            view = RunSummaryView(ds)
            df = view.summary().to_df()
            assert len(df) == 2  # run_01 and run_02


class TestQualityControlView:
    """Test QualityControlView."""

    def test_metrics_returns_cached_result(self, dataset_dir):
        from qpx.views.api import QualityControlView

        with Dataset(dataset_dir) as ds:
            view = QualityControlView(ds)
            result = view.metrics()
            assert isinstance(result, CachedQueryResult)

    def test_metrics_has_expected_columns(self, dataset_dir):
        from qpx.views.api import QualityControlView

        with Dataset(dataset_dir) as ds:
            view = QualityControlView(ds)
            df = view.metrics().to_df()
            assert "unique_peptides" in df.columns
            assert "unique_proteins" in df.columns
            assert "n_runs" in df.columns
            assert "total_features" in df.columns

    def test_metrics_counts_non_decoy_only(self, dataset_dir):
        from qpx.views.api import QualityControlView

        with Dataset(dataset_dir) as ds:
            view = QualityControlView(ds)
            df = view.metrics().to_df()
            # All 3 test features have is_decoy=False
            assert df.iloc[0]["total_features"] == 3


class TestDatasetViewProperties:
    """Test Dataset convenience properties for views."""

    def test_protein_view_property(self, dataset_dir):
        with Dataset(dataset_dir) as ds:
            assert ds.protein_view is not None
            assert hasattr(ds.protein_view, "intensity")

    def test_peptide_view_property(self, dataset_dir):
        with Dataset(dataset_dir) as ds:
            assert ds.peptide_view is not None
            assert hasattr(ds.peptide_view, "intensity")

    def test_identification_summary_property(self, dataset_dir):
        with Dataset(dataset_dir) as ds:
            assert ds.identification_summary is not None
            assert hasattr(ds.identification_summary, "summary")

    def test_run_summary_property(self, dataset_dir):
        with Dataset(dataset_dir) as ds:
            assert ds.run_summary is not None
            assert hasattr(ds.run_summary, "summary")

    def test_qc_view_property(self, dataset_dir):
        with Dataset(dataset_dir) as ds:
            assert ds.qc_view is not None
            assert hasattr(ds.qc_view, "metrics")

    def test_view_property_produces_results(self, dataset_dir):
        with Dataset(dataset_dir) as ds:
            df = ds.run_summary.summary().to_df()
            assert len(df) > 0
