"""Shared fixtures for QPX new API tests.

Provides temp directories, minimal valid Parquet files for each structure type,
and a complete dataset directory with feature, pg, sample, and run parquets.
"""

from pathlib import Path

import pytest

from qpx.writers import (
    DatasetWriter,
    FeatureWriter,
    OntologyWriter,
    PgWriter,
    ProvenanceWriter,
    PsmWriter,
    RunWriter,
    SampleWriter,
)

# Stable root paths for use by tests in subdirectories
TESTS_ROOT = Path(__file__).parent
PROJECT_ROOT = TESTS_ROOT.parent

# ---------------------------------------------------------------------------
# Test record factories -- produce minimal valid dicts for each schema
# ---------------------------------------------------------------------------


def make_feature_record(
    sequence="PEPTIDEK",
    peptidoform="PEPTIDEK",
    charge=2,
    is_decoy=False,
    calculated_mz=450.25,
    observed_mz=450.26,
    run_file_name="run_01",
    anchor_protein="P12345",
    scan=None,
    intensities=None,
):
    """Create a minimal valid feature record dict."""
    if scan is None:
        scan = [1001]
    if intensities is None:
        intensities = [{"label": "TMT126", "intensity": 1000.0}]
    return {
        "sequence": sequence,
        "peptidoform": peptidoform,
        "modifications": None,
        "charge": charge,
        "posterior_error_probability": 0.001,
        "is_decoy": is_decoy,
        "calculated_mz": calculated_mz,
        "observed_mz": observed_mz,
        "additional_scores": None,
        "predicted_rt": None,
        "run_file_name": run_file_name,
        "cv_params": None,
        "scan": scan,
        "rt": 120.5,
        "ion_mobility": None,
        "intensities": intensities,
        "additional_intensities": None,
        "pg_accessions": [{"accession": "P12345", "start": None, "end": None}],
        "anchor_protein": anchor_protein,
        "unique": True,
        "pg_global_qvalue": None,
        "pg_positions": None,
        "ion_mobility_start": None,
        "ion_mobility_stop": None,
        "gg_accessions": None,
        "gg_names": None,
        "id_run_file_name": None,
        "rt_start": None,
        "rt_stop": None,
    }


def make_psm_record(
    sequence="PEPTIDEK",
    peptidoform="PEPTIDEK",
    charge=2,
    is_decoy=False,
    calculated_mz=450.25,
    observed_mz=450.26,
    run_file_name="run_01",
    scan=None,
):
    """Create a minimal valid PSM record dict."""
    if scan is None:
        scan = [1001]
    return {
        "sequence": sequence,
        "peptidoform": peptidoform,
        "modifications": None,
        "charge": charge,
        "posterior_error_probability": 0.001,
        "is_decoy": is_decoy,
        "calculated_mz": calculated_mz,
        "observed_mz": observed_mz,
        "additional_scores": None,
        "predicted_rt": None,
        "run_file_name": run_file_name,
        "cv_params": None,
        "scan": scan,
        "rt": 120.5,
        "ion_mobility": None,
        "protein_accessions": ["P12345"],
        "cross_links": None,
        "mz_array": None,
        "intensity_array": None,
        "charge_array": None,
        "ion_type_array": None,
        "ion_mobility_array": None,
    }


def make_pg_record(
    anchor_protein="P12345",
    run_file_name="run_01",
    is_decoy=False,
    intensities=None,
):
    """Create a minimal valid PG record dict."""
    if intensities is None:
        intensities = [{"label": "TMT126", "intensity": 5000.0}]
    return {
        "pg_accessions": ["P12345", "P12346"],
        "pg_names": None,
        "gg_accessions": None,
        "gg_names": ["GENE1"],
        "anchor_protein": anchor_protein,
        "run_file_name": run_file_name,
        "global_qvalue": 0.005,
        "pg_qvalue": None,
        "intensities": intensities,
        "additional_intensities": None,
        "is_decoy": is_decoy,
        "contaminant": False,
        "peptides": [{"protein_name": "P12345", "peptide_count": 5}],
        "peptide_counts": {"unique_sequences": 3, "total_sequences": 5},
        "feature_counts": None,
        "sequence_coverage": None,
        "molecular_weight": None,
        "additional_scores": None,
        "cv_params": None,
    }


def make_sample_record(sample_accession="SAMPLE_01"):
    """Create a minimal valid sample record dict."""
    return {
        "sample_accession": sample_accession,
        "organism": "Homo sapiens",
        "organism_part": "liver",
        "disease": None,
        "cell_line": None,
        "cell_type": None,
        "sex": None,
        "age": None,
        "developmental_stage": None,
        "ancestry": None,
        "individual": None,
        "sample_description": None,
    }


def make_run_record(run_accession="assay_01", run_file_name="run_01"):
    """Create a minimal valid run record dict."""
    return {
        "run_accession": run_accession,
        "run_file_name": run_file_name,
        "samples": [
            {
                "sample_accession": "SAMPLE_01",
                "label": "TMT126",
                "biological_replicate": 1,
                "technical_replicate": 1,
            }
        ],
        "fraction": None,
        "instrument": None,
        "enzymes": None,
        "dissociation_method": None,
        "modification_parameters": None,
    }


def make_dataset_record(project_accession="PXD000001"):
    """Create a minimal valid dataset metadata record dict."""
    return {
        "project_accession": project_accession,
        "project_title": "Test Project",
        "project_description": "A test project for QPX",
        "pubmed_id": None,
        "software_name": "quantms",
        "software_version": "1.0.0",
        "creation_date": "2024-01-01T00:00:00",
        "file_checksums": None,
        "file_row_counts": None,
        "file_sizes_bytes": None,
        "total_structures": None,
        "packaged_at": None,
    }


def make_ontology_record(field_name="sequence", view="feature"):
    """Create a minimal valid ontology record dict."""
    return {
        "field_name": field_name,
        "ontology_name": "peptide sequence",
        "ontology_accession": "MS:1000888",
        "ontology_source": "MS",
        "ontology_version": "4.1.235",
        "view": view,
        "description": "The amino acid sequence of the peptide",
    }


def make_provenance_record(step_order=1):
    """Create a minimal valid provenance record dict."""
    return {
        "step_order": step_order,
        "step_category": "database_search",
        "step_name": "sequence_search",
        "tool_name": "MSFragger",
        "tool_version": "3.5",
        "tool_uri": None,
        "parameters": [{"key": "precursor_tol", "value": "20ppm"}],
        "config": None,
        "output_views": ["psm"],
    }


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def feature_parquet(tmp_path):
    """Write a minimal feature parquet with 3 rows."""
    path = tmp_path / "test.feature.parquet"
    records = [
        make_feature_record(sequence="PEPTIDEK", anchor_protein="P12345", run_file_name="run_01"),
        make_feature_record(
            sequence="ANOTHERK",
            anchor_protein="P12345",
            run_file_name="run_01",
            charge=3,
            calculated_mz=500.30,
            observed_mz=500.31,
            scan=[1002],
        ),
        make_feature_record(
            sequence="PEPTIDEK",
            anchor_protein="P67890",
            run_file_name="run_02",
            is_decoy=True,
            scan=[2001],
            intensities=[{"label": "TMT127", "intensity": 2000.0}],
        ),
    ]
    with FeatureWriter(path) as w:
        w.write_batch(records)
    return path


@pytest.fixture
def psm_parquet(tmp_path):
    """Write a minimal PSM parquet with 3 rows."""
    path = tmp_path / "test.psm.parquet"
    records = [
        make_psm_record(sequence="PEPTIDEK", is_decoy=False, run_file_name="run_01"),
        make_psm_record(
            sequence="ANOTHERK",
            is_decoy=False,
            run_file_name="run_01",
            charge=3,
            scan=[1002],
        ),
        make_psm_record(
            sequence="DECOYPEP",
            is_decoy=True,
            run_file_name="run_02",
            scan=[2001],
        ),
    ]
    with PsmWriter(path) as w:
        w.write_batch(records)
    return path


@pytest.fixture
def pg_parquet(tmp_path):
    """Write a minimal PG parquet with 3 rows."""
    path = tmp_path / "test.pg.parquet"
    records = [
        make_pg_record(anchor_protein="P12345", run_file_name="run_01"),
        make_pg_record(anchor_protein="P67890", run_file_name="run_01"),
        make_pg_record(
            anchor_protein="DECOY_P99",
            run_file_name="run_02",
            is_decoy=True,
            intensities=[{"label": "TMT127", "intensity": 3000.0}],
        ),
    ]
    with PgWriter(path) as w:
        w.write_batch(records)
    return path


@pytest.fixture
def sample_parquet(tmp_path):
    """Write a minimal sample parquet with 2 rows."""
    path = tmp_path / "test.sample.parquet"
    records = [
        make_sample_record(sample_accession="SAMPLE_01"),
        make_sample_record(sample_accession="SAMPLE_02"),
    ]
    with SampleWriter(path) as w:
        w.write_batch(records)
    return path


@pytest.fixture
def run_parquet(tmp_path):
    """Write a minimal run parquet with 2 rows."""
    path = tmp_path / "test.run.parquet"
    records = [
        make_run_record(run_accession="assay_01", run_file_name="run_01"),
        make_run_record(run_accession="assay_02", run_file_name="run_02"),
    ]
    with RunWriter(path) as w:
        w.write_batch(records)
    return path


@pytest.fixture
def dataset_parquet(tmp_path):
    """Write a minimal dataset metadata parquet with 1 row."""
    path = tmp_path / "test.dataset.parquet"
    records = [make_dataset_record()]
    with DatasetWriter(path) as w:
        w.write_batch(records)
    return path


@pytest.fixture
def ontology_parquet(tmp_path):
    """Write a minimal ontology parquet with 2 rows."""
    path = tmp_path / "test.ontology.parquet"
    records = [
        make_ontology_record(field_name="sequence", view="feature"),
        make_ontology_record(field_name="anchor_protein", view="pg"),
    ]
    with OntologyWriter(path) as w:
        w.write_batch(records)
    return path


@pytest.fixture
def provenance_parquet(tmp_path):
    """Write a minimal provenance parquet with 2 rows."""
    path = tmp_path / "test.provenance.parquet"
    records = [
        make_provenance_record(step_order=1),
        make_provenance_record(step_order=2),
    ]
    with ProvenanceWriter(path) as w:
        w.write_batch(records)
    return path


@pytest.fixture
def dataset_dir(tmp_path):
    """Create a complete dataset directory with feature, pg, sample, run, and dataset parquets.

    This uses the same run_file_name and sample relationships so that
    cross-structure queries (joins, views) work correctly.
    """
    ds_dir = tmp_path / "test_dataset"
    ds_dir.mkdir()

    # Feature file
    feature_records = [
        make_feature_record(sequence="PEPTIDEK", anchor_protein="P12345", run_file_name="run_01"),
        make_feature_record(
            sequence="ANOTHERK",
            anchor_protein="P12345",
            run_file_name="run_01",
            charge=3,
            calculated_mz=500.30,
            observed_mz=500.31,
            scan=[1002],
        ),
        make_feature_record(
            sequence="THIRDPEP",
            anchor_protein="P67890",
            run_file_name="run_02",
            scan=[2001],
            intensities=[{"label": "TMT127", "intensity": 2000.0}],
        ),
    ]
    with FeatureWriter(ds_dir / "exp.feature.parquet") as w:
        w.write_batch(feature_records)

    # PG file
    pg_records = [
        make_pg_record(anchor_protein="P12345", run_file_name="run_01"),
        make_pg_record(
            anchor_protein="P67890",
            run_file_name="run_02",
            intensities=[{"label": "TMT127", "intensity": 3000.0}],
        ),
    ]
    with PgWriter(ds_dir / "exp.pg.parquet") as w:
        w.write_batch(pg_records)

    # Sample file
    sample_records = [
        make_sample_record(sample_accession="SAMPLE_01"),
        make_sample_record(sample_accession="SAMPLE_02"),
    ]
    with SampleWriter(ds_dir / "exp.sample.parquet") as w:
        w.write_batch(sample_records)

    # Run file -- run_01 maps to SAMPLE_01/TMT126, run_02 maps to SAMPLE_02/TMT127
    run_records = [
        {
            "run_accession": "assay_01",
            "run_file_name": "run_01",
            "samples": [
                {
                    "sample_accession": "SAMPLE_01",
                    "label": "TMT126",
                    "biological_replicate": 1,
                    "technical_replicate": 1,
                }
            ],
            "fraction": None,
            "instrument": None,
            "enzymes": None,
            "dissociation_method": None,
            "modification_parameters": None,
        },
        {
            "run_accession": "assay_02",
            "run_file_name": "run_02",
            "samples": [
                {
                    "sample_accession": "SAMPLE_02",
                    "label": "TMT127",
                    "biological_replicate": 2,
                    "technical_replicate": 1,
                }
            ],
            "fraction": None,
            "instrument": None,
            "enzymes": None,
            "dissociation_method": None,
            "modification_parameters": None,
        },
    ]
    with RunWriter(ds_dir / "exp.run.parquet") as w:
        w.write_batch(run_records)

    # Dataset metadata
    with DatasetWriter(ds_dir / "exp.dataset.parquet") as w:
        w.write_batch([make_dataset_record()])

    return ds_dir
