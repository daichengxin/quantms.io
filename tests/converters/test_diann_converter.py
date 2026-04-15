"""DIA-NN converter integration tests using small example dataset.

Runs the full DIA-NN conversion pipeline once (module scope) and validates
that feature.parquet, pg.parquet, and ontology.parquet are written with
correct schemas and plausible values.
"""

import re
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "examples" / "diann" / "small"

_REPORT = EXAMPLES_DIR / "diann_report.tsv"
_SDRF = EXAMPLES_DIR / "PXD019909-DIA.sdrf.tsv"
_PG_MATRIX = EXAMPLES_DIR / "diann_report.pg_matrix.tsv"
_MZML_TSV_DIR = EXAMPLES_DIR / "mzml"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _parse_scan_number(spectrum_id: str) -> int:
    """Extract numeric scan number from a SpectrumID string like
    'controllerType=0 controllerNumber=1 scan=42'.
    """
    m = re.search(r"scan=(\d+)", str(spectrum_id))
    return int(m.group(1)) if m else 0


def _prepare_ms_info_parquet(source_dir: Path, dest_dir: Path) -> None:
    """Convert *_mzml_info.tsv files into *_ms_info.parquet files.

    The feature adapter expects ``*_ms_info.parquet`` files with columns
    ``rt``, ``scan``, and ``precursor_mz``.  The test data ships as TSV
    with different column names, so this helper bridges the gap.
    """
    for tsv_path in source_dir.glob("*_mzml_info.tsv"):
        df = pd.read_csv(tsv_path, sep="\t")
        out = pd.DataFrame(
            {
                "rt": df["Retention_Time"].astype(float),
                "scan": df["SpectrumID"].apply(_parse_scan_number).astype(int),
                "precursor_mz": pd.to_numeric(df["Exp_Mass_To_Charge"], errors="coerce"),
            }
        )
        stem = tsv_path.stem.replace("_mzml_info", "")
        parquet_name = f"{stem}_ms_info.parquet"
        pq.write_table(pa.Table.from_pandas(out), str(dest_dir / parquet_name))


# ---------------------------------------------------------------------------
# Module-scoped fixture: run conversion once, share outputs across tests
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def converted_output(tmp_path_factory):
    """Run DIA-NN conversion once for all tests in this module.

    Returns the output directory containing the generated Parquet files.
    If the test data is missing, all tests in this module are skipped.
    """
    from qpx.converters.diann.converter import DiaNNConverter

    if not _REPORT.exists():
        pytest.skip(f"Test data not found: {_REPORT}")

    output = tmp_path_factory.mktemp("diann_output")

    # Prepare ms_info parquet files from the shipped TSVs.
    mzml_parquet_dir = tmp_path_factory.mktemp("mzml_info")
    _prepare_ms_info_parquet(_MZML_TSV_DIR, mzml_parquet_dir)

    converter = DiaNNConverter(
        report_path=str(_REPORT),
        sdrf_path=str(_SDRF),
    )
    converter.convert_features(
        mzml_info_folder=str(mzml_parquet_dir),
        qvalue_threshold=0.05,
        output_folder=str(output),
        output_prefix="diann_test",
    )
    converter.convert_pg(
        pg_matrix_path=str(_PG_MATRIX),
        output_folder=str(output),
        output_prefix="diann_test",
    )
    converter.write_ontology(
        output_folder=output,
        prefix="diann_test",
    )

    return output


@pytest.fixture(scope="module")
def feature_table(converted_output):
    """Read the feature.parquet produced by the converter."""
    path = converted_output / "diann_test.feature.parquet"
    if not path.exists():
        pytest.skip("feature.parquet was not produced")
    return pq.read_table(str(path))


@pytest.fixture(scope="module")
def pg_table(converted_output):
    """Read the pg.parquet produced by the converter."""
    path = converted_output / "diann_test.pg.parquet"
    if not path.exists():
        pytest.skip("pg.parquet was not produced")
    return pq.read_table(str(path))


# ---------------------------------------------------------------------------
# Feature conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestDiaNNFeatureConversion:
    """Validate feature.parquet output from DIA-NN conversion."""

    def test_file_exists(self, converted_output):
        path = converted_output / "diann_test.feature.parquet"
        assert path.exists(), "feature.parquet was not created"

    def test_has_rows(self, feature_table):
        assert feature_table.num_rows > 0, "feature.parquet is empty"

    def test_key_columns_present(self, feature_table):
        column_names = set(feature_table.column_names)
        expected = {"sequence", "charge", "intensities", "run_file_name"}
        missing = expected - column_names
        assert not missing, f"Missing columns: {missing}"

    def test_sequence_values_are_nonempty(self, feature_table):
        sequences = feature_table.column("sequence").to_pylist()
        for seq in sequences:
            assert isinstance(seq, str) and len(seq) > 0, f"Expected non-empty string, got {seq!r}"

    def test_charge_values_are_valid(self, feature_table):
        charges = feature_table.column("charge").to_pylist()
        for charge in charges:
            assert 1 <= charge <= 10, f"Charge out of range: {charge}"

    def test_intensities_are_nonnegative(self, feature_table):
        rows = feature_table.column("intensities").to_pylist()
        for row_intensities in rows:
            assert row_intensities is not None and len(row_intensities) > 0
            for entry in row_intensities:
                assert entry["intensity"] >= 0, f"Negative intensity: {entry['intensity']}"

    def test_run_file_names_are_nonempty(self, feature_table):
        run_names = feature_table.column("run_file_name").to_pylist()
        for name in run_names:
            assert isinstance(name, str) and len(name) > 0

    def test_schema_validation(self, feature_table):
        """Verify that the output conforms to the FeatureSchema."""
        from qpx.core.data import FeatureSchema

        errors = FeatureSchema.validate(feature_table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# Protein group conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestDiaNNPgConversion:
    """Validate pg.parquet output from DIA-NN conversion."""

    def test_file_exists(self, converted_output):
        path = converted_output / "diann_test.pg.parquet"
        assert path.exists(), "pg.parquet was not created"

    def test_has_rows(self, pg_table):
        assert pg_table.num_rows > 0, "pg.parquet is empty"

    def test_key_columns_present(self, pg_table):
        column_names = set(pg_table.column_names)
        expected = {"pg_accessions", "anchor_protein", "run_file_name", "intensities"}
        missing = expected - column_names
        assert not missing, f"Missing columns: {missing}"

    def test_protein_accessions_are_nonempty(self, pg_table):
        accessions_col = pg_table.column("pg_accessions").to_pylist()
        for accessions in accessions_col:
            assert accessions is not None and len(accessions) > 0
            for acc in accessions:
                assert isinstance(acc, str) and len(acc) > 0, f"Expected non-empty accession, got {acc!r}"

    def test_anchor_protein_is_nonempty(self, pg_table):
        anchors = pg_table.column("anchor_protein").to_pylist()
        for anchor in anchors:
            assert isinstance(anchor, str) and len(anchor) > 0

    def test_run_file_names_are_nonempty(self, pg_table):
        run_names = pg_table.column("run_file_name").to_pylist()
        for name in run_names:
            assert isinstance(name, str) and len(name) > 0

    def test_intensities_structure(self, pg_table):
        rows = pg_table.column("intensities").to_pylist()
        for row_intensities in rows:
            assert row_intensities is not None and len(row_intensities) > 0
            for entry in row_intensities:
                assert "label" in entry
                assert "intensity" in entry

    def test_schema_validation(self, pg_table):
        """Verify that the output conforms to the PgSchema."""
        from qpx.core.data import PgSchema

        errors = PgSchema.validate(pg_table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# Ontology conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestDiaNNOntologyConversion:
    """Validate ontology.parquet output from DIA-NN conversion.

    The ontology file is only written when scores are discovered during
    conversion. It may legitimately not exist if no scores were tracked.
    """

    def test_ontology_file_exists_or_no_scores(self, converted_output):
        path = converted_output / "diann_test.ontology.parquet"
        if not path.exists():
            pytest.skip("ontology.parquet not written (no scores discovered)")

        table = pq.read_table(str(path))
        assert table.num_rows > 0, "ontology.parquet exists but is empty"

    def test_ontology_columns(self, converted_output):
        path = converted_output / "diann_test.ontology.parquet"
        if not path.exists():
            pytest.skip("ontology.parquet not written (no scores discovered)")

        table = pq.read_table(str(path))
        column_names = set(table.column_names)
        expected = {"field_name", "view"}
        missing = expected - column_names
        assert not missing, f"Missing ontology columns: {missing}"

    def test_schema_validation(self, converted_output):
        from qpx.core.data import OntologySchema

        path = converted_output / "diann_test.ontology.parquet"
        if not path.exists():
            pytest.skip("ontology.parquet not written (no scores discovered)")
        table = pq.read_table(str(path))
        errors = OntologySchema.validate(table)
        assert not errors, f"Schema validation errors: {errors}"
