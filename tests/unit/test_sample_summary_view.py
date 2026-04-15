"""Tests for SampleSummaryView."""

import importlib.util

import pytest

from qpx import Dataset

_has_plotly = importlib.util.find_spec("plotly") is not None
requires_plotly = pytest.mark.skipif(not _has_plotly, reason="plotly not installed")


def test_sample_summary_counts(dataset_dir):
    """peptides_per_run, peptides_per_sample, proteins_per_sample all return correct data."""
    ds = Dataset(dataset_dir)

    # peptides_per_run
    df = ds.sample_summary.peptides_per_run().to_df()
    assert "run_file_name" in df.columns
    assert "n_peptides" in df.columns
    assert len(df) == 2
    run01 = df[df["run_file_name"] == "run_01"]
    assert run01["n_peptides"].iloc[0] == 2
    run02 = df[df["run_file_name"] == "run_02"]
    assert run02["n_peptides"].iloc[0] == 1

    # peptides_per_sample
    df = ds.sample_summary.peptides_per_sample().to_df()
    assert "sample_accession" in df.columns
    assert "n_peptides" in df.columns
    assert len(df) == 2

    # proteins_per_sample
    df = ds.sample_summary.proteins_per_sample().to_df()
    assert "sample_accession" in df.columns
    assert "n_proteins" in df.columns
    assert len(df) == 2


@requires_plotly
def test_sample_summary_plots(dataset_dir):
    """Plot methods return figures; unknown metric raises ValueError."""
    ds = Dataset(dataset_dir)

    for metric in ("peptides_per_run", "peptides_per_sample", "proteins_per_sample"):
        fig = ds.sample_summary.plot(metric)
        assert fig is not None
        assert hasattr(fig, "data")

    with pytest.raises(ValueError, match="Unknown metric"):
        ds.sample_summary.plot("invalid_metric")
