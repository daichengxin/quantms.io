"""MuData export — assemble QPX dataset into a MuData object.

Each modality (precursors, proteins, expression, differential) becomes
an AnnData layer inside the MuData container.  Feature-to-protein
mapping is stored as a boolean sparse adjacency matrix in
``mdata.varp["feature_mapping"]``.
"""

from __future__ import annotations

import logging
from pathlib import Path, PurePosixPath
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy import sparse

if TYPE_CHECKING:
    from qpx.core.engine import DuckDBEngine
    from qpx.dataset import Dataset

logger = logging.getLogger(__name__)

_VALID_MODALITIES = {"precursors", "proteins", "expression", "differential"}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _pivot_to_sparse(
    df: pd.DataFrame,
    row_col: str,
    col_col: str,
    value_col: str,
    row_index: pd.Index,
    col_index: pd.Index,
) -> sparse.csr_matrix:
    """Pivot a long-form DataFrame into a sparse CSR matrix.

    Parameters
    ----------
    df : DataFrame
        Long-form data with at least *row_col*, *col_col*, *value_col*.
    row_col, col_col, value_col : str
        Column names for rows, columns, and values.
    row_index, col_index : pd.Index
        Pre-built indices that define the matrix shape.

    Returns
    -------
    scipy.sparse.csr_matrix
    """
    row_pos = row_index.get_indexer(df[row_col])
    col_pos = col_index.get_indexer(df[col_col])

    # Drop entries that did not match (-1)
    mask = (row_pos >= 0) & (col_pos >= 0)
    data = df[value_col].values[mask].astype(np.float32)
    row_pos = row_pos[mask]
    col_pos = col_pos[mask]

    mat = sparse.coo_matrix(
        (data, (row_pos, col_pos)),
        shape=(len(row_index), len(col_index)),
    )
    return mat.tocsr()


_LABEL_FIELD_QUERIES: dict[tuple[str, str], str] = {
    ("channel", "feature"): "SELECT i.channel FROM feature, UNNEST(intensities) AS _t(i) LIMIT 1",
    ("channel", "pg"): "SELECT i.channel FROM pg, UNNEST(intensities) AS _t(i) LIMIT 1",
    ("label", "feature"): "SELECT i.label FROM feature, UNNEST(intensities) AS _t(i) LIMIT 1",
    ("label", "pg"): "SELECT i.label FROM pg, UNNEST(intensities) AS _t(i) LIMIT 1",
}


_PRECURSOR_QUERIES: dict[str, str] = {
    "channel": """
    SELECT f.run_file_name,
           f.peptidoform || '|' || CAST(f.charge AS VARCHAR) AS precursor_id,
           i.intensity
    FROM feature f, UNNEST(f.intensities) AS _t(i)
    WHERE i.channel = $1
    """,
    "label": """
    SELECT f.run_file_name,
           f.peptidoform || '|' || CAST(f.charge AS VARCHAR) AS precursor_id,
           i.intensity
    FROM feature f, UNNEST(f.intensities) AS _t(i)
    WHERE i.label = $1
    """,
}

_PROTEIN_QUERIES: dict[str, str] = {
    "channel": """
    SELECT pg.run_file_name,
           pg.anchor_protein,
           i.intensity
    FROM pg, UNNEST(pg.intensities) AS _t(i)
    WHERE i.channel = $1
    """,
    "label": """
    SELECT pg.run_file_name,
           pg.anchor_protein,
           i.intensity
    FROM pg, UNNEST(pg.intensities) AS _t(i)
    WHERE i.label = $1
    """,
}


_TYPEOF_QUERIES: dict[str, str] = {
    "feature": "SELECT typeof(intensities) FROM feature LIMIT 1",
    "pg": "SELECT typeof(intensities) FROM pg LIMIT 1",
}


def _detect_label_field(engine: DuckDBEngine, table: str = "feature") -> str:
    """Return ``"channel"`` or ``"label"`` by inspecting the intensities struct type."""
    try:
        row = engine.execute(_TYPEOF_QUERIES[table]).fetchone()
        if row:
            type_str = row[0].lower()
            if "channel" in type_str and "label" not in type_str:
                return "channel"
        return "label"
    except Exception:
        logger.debug("Could not detect label field from %s; defaulting to 'label'", table)
        return "label"


def _detect_intensity_label(engine: DuckDBEngine, table: str = "feature") -> str:
    """Auto-detect the first intensity label from the given table."""
    label_field = _detect_label_field(engine, table)

    query = _LABEL_FIELD_QUERIES.get((label_field, table))
    if query is None:
        raise ValueError(f"Invalid label_field={label_field!r} or table={table!r}")
    try:
        row = engine.execute(query).fetchone()
    except Exception as exc:
        raise ValueError(f"Failed to read intensity labels from {table}: {exc}") from exc
    if row is None:
        raise ValueError(f"No intensity labels found in {table}")
    return row[0]


def _load_run_obs(engine: DuckDBEngine, run_index: pd.Index) -> pd.DataFrame:
    """Load scalar run metadata as an obs DataFrame.

    Returns a DataFrame indexed by *run_index* with columns from the
    run table (run_accession, etc.) and flattened sample info.
    """
    try:
        df = engine.execute(
            """
            SELECT r.run_file_name,
                   r.run_accession,
                   rs.sample_accession,
                   rs.label,
                   rs.biological_replicate,
                   rs.technical_replicate,
                   r.fraction,
                   r.instrument
            FROM run r, UNNEST(r.samples) AS _t(rs)
            """
        ).fetchdf()
    except Exception:
        # If run table not available, return empty obs
        return pd.DataFrame(index=run_index)

    if df.empty:
        return pd.DataFrame(index=run_index)

    # Keep first row per run_file_name (one sample per run for obs)
    df = df.drop_duplicates(subset="run_file_name", keep="first")
    df = df.set_index("run_file_name")

    # Reindex to match run_index; if direct match fails (e.g. extension
    # mismatch between run.parquet and feature.parquet), try stem-based
    # matching so that "BSA1_F1" matches "BSA1_F1.mzML".
    df = df.reindex(run_index)
    if df.isna().all(axis=None) and not df.empty:
        df_stems = df.index.map(lambda n: PurePosixPath(n).stem)
        if not df_stems.duplicated().any():
            # Re-query with original index, map via stems
            df2 = engine.execute(
                """
                SELECT r.run_file_name,
                       r.run_accession,
                       rs.sample_accession,
                       rs.label,
                       rs.biological_replicate,
                       rs.technical_replicate,
                       r.fraction,
                       r.instrument
                FROM run r, UNNEST(r.samples) AS _t(rs)
                """
            ).fetchdf()
            df2 = df2.drop_duplicates(subset="run_file_name", keep="first")
            df2["_stem"] = df2["run_file_name"].map(lambda n: PurePosixPath(n).stem)
            df2 = df2.set_index("_stem").drop(columns=["run_file_name"])
            df2.index.name = "run_file_name"
            # Map stems back to original run_index names
            df2 = df2.reindex([PurePosixPath(n).stem for n in run_index])
            df2.index = run_index
            df = df2

    # Fill remaining NaN in string columns to avoid h5py serialization errors
    str_cols = df.select_dtypes(include=["object"]).columns
    df[str_cols] = df[str_cols].fillna("")
    return df


def _attach_uns_metadata(engine: DuckDBEngine, mdata) -> None:
    """Attach project metadata from dataset.parquet to MuData.uns."""
    try:
        df = engine.execute("SELECT * FROM dataset LIMIT 1").fetchdf()
    except Exception:
        logger.debug("No dataset table found; skipping uns metadata.")
        return

    if df.empty:
        return

    row = df.iloc[0]
    for col in df.columns:
        val = row[col]
        if val is None:
            continue
        # pd.isna() catches np.nan, pd.NA, pd.NaT on scalars; guarded by is_scalar
        # so that dicts/lists/arrays do not raise or get misclassified.
        if pd.api.types.is_scalar(val) and pd.isna(val):
            continue
        mdata.uns[col] = val


# ---------------------------------------------------------------------------
# Modality builders
# ---------------------------------------------------------------------------


def _build_precursor_adata(engine: DuckDBEngine, intensity_label: str, label_field: str = "label"):
    """Build AnnData from feature.parquet.

    obs = runs, var = peptidoform|charge, X = sparse intensity matrix.
    """
    import anndata as ad

    df = engine.execute(_PRECURSOR_QUERIES[label_field], [intensity_label]).fetchdf()

    if df.empty:
        return ad.AnnData()

    runs = pd.Index(sorted(df["run_file_name"].unique()), name="run_file_name")
    precursors = pd.Index(sorted(df["precursor_id"].unique()), name="precursor_id")

    X = _pivot_to_sparse(df, "run_file_name", "precursor_id", "intensity", runs, precursors)

    obs = _load_run_obs(engine, runs)
    var = pd.DataFrame(index=precursors)

    return ad.AnnData(X=X, obs=obs, var=var)


def _build_protein_adata(engine: DuckDBEngine, intensity_label: str, label_field: str = "label"):
    """Build AnnData from pg.parquet.

    obs = runs, var = anchor_protein, X = sparse intensity matrix.
    var includes gene_name column.
    """
    import anndata as ad

    df = engine.execute(_PROTEIN_QUERIES[label_field], [intensity_label]).fetchdf()

    if df.empty:
        return ad.AnnData()

    runs = pd.Index(sorted(df["run_file_name"].unique()), name="run_file_name")
    proteins = pd.Index(sorted(df["anchor_protein"].unique()), name="anchor_protein")

    X = _pivot_to_sparse(df, "run_file_name", "anchor_protein", "intensity", runs, proteins)

    obs = _load_run_obs(engine, runs)

    # Gene names
    gene_sql = """
    SELECT DISTINCT pg.anchor_protein,
           pg.gg_names[1] AS gene_name
    FROM pg
    """
    gene_df = engine.execute(gene_sql).fetchdf()
    gene_df = gene_df.drop_duplicates(subset="anchor_protein", keep="first")
    gene_df = gene_df.set_index("anchor_protein").reindex(proteins)

    var = pd.DataFrame({"gene_name": gene_df["gene_name"].values}, index=proteins)

    return ad.AnnData(X=X, obs=obs, var=var)


def _build_feature_mapping(engine: DuckDBEngine, mdata) -> sparse.csr_matrix:
    """Build boolean sparse adjacency linking precursor IDs to protein IDs.

    Returns a symmetric CSR matrix of shape (N, N) where N is the total
    number of var_names across all modalities in mdata.  Non-zero entries
    indicate a precursor-protein relationship.
    """
    if "precursors" not in mdata.mod or "proteins" not in mdata.mod:
        raise ValueError("Both 'precursors' and 'proteins' modalities are required for feature mapping.")

    precursor_names = mdata.mod["precursors"].var_names
    protein_names = mdata.mod["proteins"].var_names

    # Query mapping from feature table
    sql = """
    SELECT DISTINCT
           f.peptidoform || '|' || CAST(f.charge AS VARCHAR) AS precursor_id,
           f.anchor_protein
    FROM feature f
    """
    df = engine.execute(sql).fetchdf()

    if df.empty:
        n = len(precursor_names) + len(protein_names)
        return sparse.csr_matrix((n, n), dtype=bool)

    # Build combined index
    all_names = pd.Index(list(precursor_names) + list(protein_names))

    precursor_pos = all_names.get_indexer(df["precursor_id"])
    protein_pos = all_names.get_indexer(df["anchor_protein"])

    mask = (precursor_pos >= 0) & (protein_pos >= 0)
    precursor_pos = precursor_pos[mask]
    protein_pos = protein_pos[mask]

    n = len(all_names)
    data = np.ones(len(precursor_pos) * 2, dtype=bool)
    rows = np.concatenate([precursor_pos, protein_pos])
    cols = np.concatenate([protein_pos, precursor_pos])

    mat = sparse.coo_matrix((data, (rows, cols)), shape=(n, n))
    return mat.tocsr()


def _build_expression_adata(ds_path: Path):
    """Auto-detect and load ``.pe.h5ad`` or ``.pe.zarr``.

    Returns None if not found.
    """
    import anndata as ad

    # Try h5ad first
    candidates = list(ds_path.glob("*.pe.h5ad"))
    if candidates:
        return ad.read_h5ad(candidates[0])

    # Try zarr
    candidates = list(ds_path.glob("*.pe.zarr"))
    if candidates:
        return ad.read_zarr(candidates[0])

    return None


def _build_differential_adata(ds_path: Path):
    """Auto-detect and load ``.de.h5ad`` or ``.de.zarr``.

    Restructures ``uns["de_results"]`` into matrix form:
    obs = contrasts, var = proteins, X = log2FC,
    layers = {pvals_adj, scores, pvals, se}.

    Returns None if not found.
    """
    import anndata as ad

    # Try h5ad first
    candidates = list(ds_path.glob("*.de.h5ad"))
    if not candidates:
        candidates = list(ds_path.glob("*.de.zarr"))
        if not candidates:
            return None
        adata = ad.read_zarr(candidates[0])
    else:
        adata = ad.read_h5ad(candidates[0])

    # If it has de_results in uns, restructure into matrix form
    if "de_results" in adata.uns:
        de = adata.uns["de_results"]
        if isinstance(de, pd.DataFrame) and not de.empty:
            return _reshape_de_results(de)

    # Already in matrix form or no de_results — return as-is
    return adata


def _reshape_de_results(de: pd.DataFrame):
    """Reshape long-form DE results into AnnData matrix.

    Expects columns: contrast, protein, log2FC and optionally
    pvals_adj, scores, pvals, se.
    """
    import anndata as ad

    # Identify contrast and protein columns
    contrast_col = "contrast" if "contrast" in de.columns else de.columns[0]
    protein_col = "protein" if "protein" in de.columns else de.columns[1]

    contrasts = pd.Index(sorted(de[contrast_col].unique()), name="contrast")
    proteins = pd.Index(sorted(de[protein_col].unique()), name="protein")

    if "log2FC" not in de.columns:
        raise ValueError(f"DE results missing required 'log2FC' column. Found columns: {list(de.columns)}")

    X = _pivot_to_sparse(de, contrast_col, protein_col, "log2FC", contrasts, proteins)

    layers = {}
    for col in ("pvals_adj", "scores", "pvals", "se"):
        if col in de.columns:
            layers[col] = _pivot_to_sparse(de, contrast_col, protein_col, col, contrasts, proteins)

    return ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=contrasts),
        var=pd.DataFrame(index=proteins),
        layers=layers,
    )


# ---------------------------------------------------------------------------
# Top-level assembler
# ---------------------------------------------------------------------------


def _try_build_modality(name: str, builder, mod: dict) -> None:
    """Call *builder*, add the result to *mod* if non-empty."""
    try:
        adata = builder()
        if adata is not None and (not hasattr(adata, "n_obs") or adata.n_obs > 0):
            mod[name] = adata
    except Exception as exc:
        logger.warning("Failed to build %s modality: %s", name, exc)


def build_mudata(
    dataset: Dataset,
    intensity_label: str | None = None,
    modalities: list[str] | None = None,
):
    """Assemble QPX dataset into a MuData object.

    Parameters
    ----------
    dataset : Dataset
        An open QPX Dataset instance.
    intensity_label : str, optional
        Intensity label to extract (e.g. "TMT126").  Auto-detected from
        feature.parquet when omitted.
    modalities : list of str, optional
        Which modalities to include.  Valid values:
        ``"precursors"``, ``"proteins"``, ``"expression"``, ``"differential"``.
        Defaults to all available.

    Returns
    -------
    mudata.MuData
    """
    import mudata as mu

    engine = dataset._engine
    ds_path = Path(dataset.path)

    if modalities is not None:
        unknown = set(modalities) - _VALID_MODALITIES
        if unknown:
            raise ValueError(f"Unknown modalities: {unknown}. Valid: {sorted(_VALID_MODALITIES)}")
        requested = set(modalities)
    else:
        requested = _VALID_MODALITIES.copy()

    feat_label_field = _detect_label_field(engine, "feature")
    feat_intensity = intensity_label or _detect_intensity_label(engine, "feature")

    pg_label_field = _detect_label_field(engine, "pg")
    try:
        pg_intensity = intensity_label or _detect_intensity_label(engine, "pg")
    except ValueError:
        pg_intensity = feat_intensity

    mod: dict = {}
    builders = {
        "precursors": lambda: _build_precursor_adata(engine, feat_intensity, feat_label_field),
        "proteins": lambda: _build_protein_adata(engine, pg_intensity, pg_label_field),
        "expression": lambda: _build_expression_adata(ds_path),
        "differential": lambda: _build_differential_adata(ds_path),
    }
    for name, builder in builders.items():
        if name in requested:
            _try_build_modality(name, builder, mod)

    mdata = mu.MuData(mod)

    # Attach feature mapping if both precursors and proteins are present
    if "precursors" in mod and "proteins" in mod:
        try:
            mapping = _build_feature_mapping(engine, mdata)
            mdata.varp["feature_mapping"] = mapping
        except Exception as exc:
            logger.warning("Failed to build feature mapping: %s", exc)

    # Attach dataset-level metadata
    _attach_uns_metadata(engine, mdata)

    # Sanitize NaN in object columns across all modalities to prevent
    # h5py serialization errors (can't write NaN in vlen-string arrays).
    for adata in mdata.mod.values():
        for frame in (adata.obs, adata.var):
            str_cols = frame.select_dtypes(include=["object"]).columns
            if len(str_cols):
                frame[str_cols] = frame[str_cols].fillna("")

    return mdata
