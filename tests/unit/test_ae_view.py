"""Tests for AbsoluteExpressionView (PCA from .ae.h5ad)."""

import pytest

plotly = pytest.importorskip("plotly", reason="plotly required for plotting tests")
ad = pytest.importorskip("anndata", reason="anndata required for AE view tests")

import numpy as np


@pytest.fixture
def ae_dataset_dir(dataset_dir):
    """Extend dataset_dir with a minimal .ae.h5ad file."""
    import anndata as ad
    import pandas as pd

    X = np.array([[100.0, 200.0, 300.0], [400.0, 500.0, 600.0]], dtype=np.float32)

    obs = pd.DataFrame(
        {
            "organism": ["Homo sapiens", "Homo sapiens"],
            "disease": ["healthy", "cancer"],
        },
        index=["SAMPLE_01", "SAMPLE_02"],
    )
    var = pd.DataFrame(index=["P12345", "P67890", "P11111"])

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.layers["ibaq"] = X * 2
    adata.layers["ibaq_log2"] = np.log2(np.where(X > 0, X, np.nan)).astype(np.float32)

    prefix = dataset_dir.name
    adata.write_h5ad(dataset_dir / f"{prefix}.ae.h5ad")
    return dataset_dir


@pytest.mark.timeout(120)
def test_pca_basic_and_options(ae_dataset_dir):
    """PCA: basic, with color, with layer, with fillna."""
    from qpx import Dataset

    ds = Dataset(ae_dataset_dir)

    # Basic PCA
    fig = ds.ae_view.pca()
    assert fig is not None
    assert len(fig.data) >= 1

    # With color
    fig = ds.ae_view.pca(color_by="disease")
    assert len(fig.data) == 2

    # With layer
    fig = ds.ae_view.pca(layer="ibaq")
    assert fig is not None

    # With fillna (ibaq_log2 has NaN values)
    fig = ds.ae_view.pca(layer="ibaq_log2", fillna="median")
    assert fig is not None


def test_pca_no_h5ad_raises(dataset_dir):
    """PCA raises FileNotFoundError when no .ae.h5ad exists."""
    from qpx import Dataset

    ds = Dataset(dataset_dir)
    with pytest.raises((FileNotFoundError, ImportError)):
        ds.ae_view.pca()
