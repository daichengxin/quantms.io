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


# ---------------------------------------------------------------------------
# PG grouped_runs field (refactor: run_file_name scalar -> grouped_runs list<string>)
# ---------------------------------------------------------------------------


def test_pg_schema_grouped_runs_is_list_string():
    """The pg schema keys on `grouped_runs` typed as list<string>."""
    from qpx.core.data.pg import PgSchema

    arrow_schema = PgSchema.get_arrow_schema()
    assert "run_file_name" not in arrow_schema.names
    field = arrow_schema.field("grouped_runs")
    assert pa.types.is_list(field.type)
    assert pa.types.is_string(field.type.value_type)
    assert tuple(PgSchema._primary_key) == ("anchor_protein", "grouped_runs", "label")


def test_pg_grouped_runs_roundtrip(pg_parquet):
    """A round-tripped pg record exposes `grouped_runs` as a single-element list."""
    pg = PG.from_file(pg_parquet)
    df = pg.to_df()
    assert "grouped_runs" in df.columns
    assert "run_file_name" not in df.columns
    row = df[df["anchor_protein"] == "P12345"].iloc[0]
    assert list(row["grouped_runs"]) == ["run_01"]


def test_pg_by_run_filters_via_grouped_runs(pg_parquet):
    """PG.by_run matches on list membership of the grouped_runs list."""
    pg = PG.from_file(pg_parquet)
    assert pg.by_run("run_01").count() == 2
    assert pg.by_run("run_02").count() == 1
    assert pg.by_run("missing").count() == 0


# ---------------------------------------------------------------------------
# grouped_runs referential invariant (pg.grouped_runs subset of run.run_file_name)
# ---------------------------------------------------------------------------


def _write_pg_run_dataset(
    ds_dir,
    pg_grouped_runs,
    run_sample_accessions=("SAMPLE_01", "SAMPLE_01"),
):
    """Write a minimal pg + run dataset; pg rows use the given grouped_runs lists."""
    from qpx.writers import DatasetWriter, PgWriter, RunWriter
    from tests.conftest import (
        make_dataset_record,
        make_pg_record,
        make_run_record,
    )

    ds_dir.mkdir(parents=True, exist_ok=True)
    pg_records = []
    for idx, grouped in enumerate(pg_grouped_runs):
        rec = make_pg_record(anchor_protein=f"P{idx:05d}")
        rec["grouped_runs"] = grouped
        pg_records.append(rec)
    with PgWriter(ds_dir / "exp.pg.parquet") as w:
        w.write_batch(pg_records)
    run_records = [
        make_run_record(run_accession="assay_01", run_file_name="run_01"),
        make_run_record(run_accession="assay_02", run_file_name="run_02"),
    ]
    for record, sample_accession in zip(
        run_records,
        run_sample_accessions,
        strict=True,
    ):
        record["samples"][0]["sample_accession"] = sample_accession
    with RunWriter(ds_dir / "exp.run.parquet") as w:
        w.write_batch(run_records)
    with DatasetWriter(ds_dir / "exp.dataset.parquet") as w:
        w.write_batch([make_dataset_record()])
    return ds_dir


def test_grouped_runs_dangling_reference_is_error(tmp_path):
    """A pg.grouped_runs value with no matching run.run_file_name is an ERROR."""
    from qpx.dataset import Dataset

    ds_dir = _write_pg_run_dataset(
        tmp_path / "bad",
        pg_grouped_runs=[["run_01"], ["run_99"]],  # run_99 does not exist
    )
    with Dataset(ds_dir) as ds:
        results = ds.validate(structures=["pg"])

    pg_result = results["pg"]
    dangling = [i for i in pg_result.issues if i.check == "dangling_grouped_run"]
    assert len(dangling) == 1
    assert dangling[0].severity == "error"
    assert "run_99" in dangling[0].message
    assert not pg_result.is_valid  # the error must fail validation


def test_grouped_runs_all_valid_passes_invariant(tmp_path):
    """When every grouped_runs value is a real run_file_name, no dangling error."""
    from qpx.dataset import Dataset

    ds_dir = _write_pg_run_dataset(
        tmp_path / "good",
        pg_grouped_runs=[["run_01"], ["run_02"], ["run_01", "run_02"]],
    )
    with Dataset(ds_dir) as ds:
        results = ds.validate(structures=["pg"])

    dangling = [i for i in results["pg"].issues if i.check == "dangling_grouped_run"]
    assert dangling == []


def test_pg_auto_join_run_expands_grouped_runs(tmp_path):
    """The documented PG-to-Run join expands every grouped run."""
    from qpx.dataset import Dataset

    ds_dir = _write_pg_run_dataset(
        tmp_path / "join",
        pg_grouped_runs=[["run_01", "run_02"]],
    )
    with Dataset(ds_dir) as ds:
        joined = ds.pg.join(ds.run).to_df()

    assert len(joined) == 2
    assert set(joined["run_file_name"]) == {"run_01", "run_02"}


def test_grouped_runs_label_resolving_to_multiple_samples_is_error(tmp_path):
    """All fractions in one PG unit must resolve a label to the same sample."""
    from qpx.dataset import Dataset

    ds_dir = _write_pg_run_dataset(
        tmp_path / "ambiguous",
        pg_grouped_runs=[["run_01", "run_02"]],
        run_sample_accessions=("SAMPLE_01", "SAMPLE_02"),
    )
    with Dataset(ds_dir) as ds:
        results = ds.validate(structures=["pg"])

    ambiguous = [issue for issue in results["pg"].issues if issue.check == "ambiguous_grouped_run_sample"]
    assert len(ambiguous) == 1
    assert ambiguous[0].severity == "error"
    assert "2 samples" in ambiguous[0].message
    assert not results["pg"].is_valid
