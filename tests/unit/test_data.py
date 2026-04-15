"""Tests for QPX data layer: BaseStructure, Feature, PSM, PG and domain methods."""

import pandas as pd
import pyarrow as pa

from qpx.core.data.dataset import DatasetMeta
from qpx.core.data.feature import Feature
from qpx.core.data.ontology import Ontology
from qpx.core.data.pg import PG
from qpx.core.data.provenance import Provenance
from qpx.core.data.psm import PSM
from qpx.core.data.run import Run
from qpx.core.data.sample import Sample

# ---------------------------------------------------------------------------
# Parameterized from_file across all structures
# ---------------------------------------------------------------------------


def test_structures_from_file(
    feature_parquet,
    psm_parquet,
    pg_parquet,
    sample_parquet,
    run_parquet,
    dataset_parquet,
    ontology_parquet,
    provenance_parquet,
):
    """All structures open from Parquet, return correct type and count."""
    cases = [
        (feature_parquet, Feature, 3),
        (psm_parquet, PSM, 3),
        (pg_parquet, PG, 3),
        (sample_parquet, Sample, 2),
        (run_parquet, Run, 2),
        (dataset_parquet, DatasetMeta, 1),
        (ontology_parquet, Ontology, 2),
        (provenance_parquet, Provenance, 2),
    ]
    for path, cls, expected_count in cases:
        obj = cls.from_file(path)
        assert isinstance(obj, cls)
        assert obj.count() == expected_count


def test_feature_filter_and_query_ops(feature_parquet):
    """filter, select, limit, count, len, and chaining all work correctly."""
    feat = Feature.from_file(feature_parquet)

    # count and len
    assert feat.count() == 3
    assert len(feat) == 3

    # filter
    filtered = feat.filter("charge = 2")
    assert isinstance(filtered, Feature)
    assert filtered.count() == 2

    # by_protein / by_run domain filters
    assert feat.by_protein("P12345").count() == 2
    assert feat.by_run("run_01").count() == 2

    # filter chaining
    assert feat.filter("charge = 2").filter("is_decoy = false").count() == 1

    # select
    selected = feat.select("sequence", "charge")
    df = selected.to_df()
    assert set(df.columns) == {"sequence", "charge"}
    assert len(df) == 3

    # limit
    assert feat.limit(2).count() == 2

    # count after filter
    assert feat.filter("is_decoy = true").count() == 1


def test_feature_materialization(feature_parquet):
    """to_df returns DataFrame, to_arrow returns Arrow Table."""
    feat = Feature.from_file(feature_parquet)
    df = feat.to_df()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 3

    table = feat.to_arrow()
    assert isinstance(table, pa.Table)
    assert table.num_rows == 3


def test_feature_domain_methods(feature_parquet):
    """unique_proteins, peptide_intensities, file_metadata, and repr."""
    feat = Feature.from_file(feature_parquet)

    # unique_proteins
    assert set(feat.unique_proteins()) == {"P12345", "P67890"}

    # peptide_intensities
    result = feat.peptide_intensities()
    df = result.to_df()
    assert "label" in df.columns
    assert "intensity" in df.columns
    assert len(df) == 3

    # file_metadata
    meta = feat.file_metadata
    assert isinstance(meta, dict)
    assert meta["file_type"] == "feature_file"

    # repr
    r = repr(feat)
    assert "Feature" in r
    assert "rows=3" in r


def test_immutability(feature_parquet):
    """filter, select, limit return new instances without mutating the original."""
    feat = Feature.from_file(feature_parquet)
    original_count = feat.count()

    filtered = feat.filter("charge > 10")
    assert filtered.count() == 0
    assert feat.count() == original_count

    selected = feat.select("sequence")
    assert selected is not feat

    limited = feat.limit(1)
    assert limited.count() == 1
    assert feat.count() == original_count


def test_iter_batches(feature_parquet):
    """iter_batches partitions correctly with different batch sizes."""
    feat = Feature.from_file(feature_parquet)

    # Large batch: one batch with all distinct runs
    batches = list(feat.iter_batches(partition_by="run_file_name", batch_size=10))
    assert len(batches) == 1
    keys, df = batches[0]
    assert set(keys) == {"run_01", "run_02"}
    assert len(df) == 3

    # Small batch: one batch per distinct value
    batches = list(feat.iter_batches(partition_by="run_file_name", batch_size=1))
    assert len(batches) == 2
    for keys, df in batches:
        assert len(keys) == 1
        assert isinstance(df, pd.DataFrame)
        run_values = df["run_file_name"].unique()
        assert set(run_values) == set(keys)


def test_order_by(feature_parquet):
    """order_by sorts ascending and descending correctly."""
    feat = Feature.from_file(feature_parquet)

    charges_asc = feat.order_by("charge").to_df()["charge"].tolist()
    assert charges_asc == sorted(charges_asc)

    charges_desc = feat.order_by("charge", desc=True).to_df()["charge"].tolist()
    assert charges_desc == sorted(charges_desc, reverse=True)


def test_feature_join_run(feature_parquet, run_parquet):
    """Feature.join(Run) joins on run_file_name."""
    feat = Feature.from_file(feature_parquet)
    feat._engine.register_parquet("run", run_parquet)
    run = Run(engine=feat._engine, table_name="run", file_path=run_parquet)

    joined = feat.join(run, on="run_file_name")
    df = joined.to_df()
    assert len(df) == 3
    assert "run_accession" in df.columns
