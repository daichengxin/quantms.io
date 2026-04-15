"""View plotting tests."""

import pytest

plotly = pytest.importorskip("plotly", reason="plotly required for plotting tests")


@pytest.fixture
def ds(dataset_dir):
    import qpx

    d = qpx.open_dataset(str(dataset_dir))
    yield d
    d.close()


class TestIdentificationSummaryPlot:
    def test_plot_returns_figure(self, ds):
        import plotly.graph_objects as go

        fig = ds.identification_summary.plot()
        assert isinstance(fig, go.Figure)


class TestRunSummaryPlot:
    def test_plot_returns_figure(self, ds):
        import plotly.graph_objects as go

        fig = ds.run_summary.plot()
        assert isinstance(fig, go.Figure)


class TestQCViewPlot:
    def test_plot_returns_figure(self, ds):
        import plotly.graph_objects as go

        fig = ds.qc_view.plot()
        assert isinstance(fig, go.Figure)


class TestPlotSave:
    def test_save_to_html(self, ds, tmp_path):
        fig = ds.identification_summary.plot()
        out = tmp_path / "plot.html"
        fig.write_html(str(out))
        assert out.exists()
        assert out.stat().st_size > 0
