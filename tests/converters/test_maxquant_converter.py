"""MaxQuant converter integration tests using small example dataset.

Runs the full MaxQuant conversion pipeline once (module scope) and validates
that psm.parquet, feature.parquet, and ontology.parquet are written with
correct schemas and plausible values.
"""

from pathlib import Path

import pyarrow.parquet as pq
import pytest

EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "examples" / "maxquant" / "maxquant_simple"

_MSMS = EXAMPLES_DIR / "msms.txt"
_EVIDENCE = EXAMPLES_DIR / "evidence.txt"
_SDRF = EXAMPLES_DIR / "sdrf.tsv"

_PREFIX = "maxquant_test"


# ---------------------------------------------------------------------------
# Module-scoped fixture: run conversion once, share outputs across tests
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def converted_output(tmp_path_factory):
    """Run MaxQuant conversion once for all tests in this module.

    Returns the output directory containing the generated Parquet files.
    If the test data is missing, all tests in this module are skipped.
    """
    from qpx.converters.maxquant.converter import MaxQuantConverter

    if not _MSMS.exists() or not _EVIDENCE.exists():
        pytest.skip(f"Test data not found: {EXAMPLES_DIR}")

    output = tmp_path_factory.mktemp("maxquant_output")

    converter = MaxQuantConverter()
    converter.convert(
        output_folder=output,
        msms_file=_MSMS,
        evidence_file=_EVIDENCE,
        sdrf_file=_SDRF if _SDRF.exists() else None,
        output_prefix=_PREFIX,
    )

    return output


@pytest.fixture(scope="module")
def psm_table(converted_output):
    """Read the psm.parquet produced by the converter."""
    path = converted_output / f"{_PREFIX}.psm.parquet"
    if not path.exists():
        pytest.skip("psm.parquet was not produced")
    return pq.read_table(str(path))


@pytest.fixture(scope="module")
def feature_table(converted_output):
    """Read the feature.parquet produced by the converter."""
    path = converted_output / f"{_PREFIX}.feature.parquet"
    if not path.exists():
        pytest.skip("feature.parquet was not produced")
    return pq.read_table(str(path))


# ---------------------------------------------------------------------------
# PSM conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestMaxQuantPsmConversion:
    """Validate psm.parquet output from MaxQuant conversion."""

    def test_file_exists(self, converted_output):
        path = converted_output / f"{_PREFIX}.psm.parquet"
        if not path.exists():
            raise AssertionError("psm.parquet was not created")

    def test_has_rows(self, psm_table):
        if psm_table.num_rows == 0:
            raise AssertionError("psm.parquet is empty")

    def test_key_columns_present(self, psm_table):
        column_names = set(psm_table.column_names)
        expected = {"sequence", "charge", "calculated_mz", "run_file_name", "scan"}
        missing = expected - column_names
        if missing:
            raise AssertionError(f"Missing columns: {missing}")

    def test_sequence_values_are_nonempty(self, psm_table):
        sequences = psm_table.column("sequence").to_pylist()
        for seq in sequences:
            if not (isinstance(seq, str) and len(seq) > 0):
                raise AssertionError(f"Expected non-empty string, got {seq!r}")

    def test_charge_values_are_valid(self, psm_table):
        charges = psm_table.column("charge").to_pylist()
        for charge in charges:
            if not 1 <= charge <= 10:
                raise AssertionError(f"Charge out of range: {charge}")

    def test_calculated_mz_values_are_nonnegative(self, psm_table):
        mz_values = psm_table.column("calculated_mz").to_pylist()
        for mz in mz_values:
            if mz is None or mz < 0:
                raise AssertionError(f"Invalid calculated_mz value: {mz}")

    def test_run_file_names_are_nonempty(self, psm_table):
        run_names = psm_table.column("run_file_name").to_pylist()
        for name in run_names:
            if not (isinstance(name, str) and len(name) > 0):
                raise AssertionError(f"Expected non-empty run name, got {name!r}")

    def test_schema_validation(self, psm_table):
        """Verify that the output conforms to the PsmSchema."""
        from qpx.core.data import PsmSchema

        errors = PsmSchema.validate(psm_table)
        if errors:
            raise AssertionError(f"Schema validation errors: {errors}")


# ---------------------------------------------------------------------------
# Feature conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestMaxQuantFeatureConversion:
    """Validate feature.parquet output from MaxQuant conversion."""

    def test_file_exists(self, converted_output):
        path = converted_output / f"{_PREFIX}.feature.parquet"
        if not path.exists():
            raise AssertionError("feature.parquet was not created")

    def test_has_rows(self, feature_table):
        if feature_table.num_rows == 0:
            raise AssertionError("feature.parquet is empty")

    def test_key_columns_present(self, feature_table):
        column_names = set(feature_table.column_names)
        expected = {
            "sequence",
            "charge",
            "calculated_mz",
            "intensities",
            "run_file_name",
        }
        missing = expected - column_names
        if missing:
            raise AssertionError(f"Missing columns: {missing}")

    def test_sequence_values_are_nonempty(self, feature_table):
        sequences = feature_table.column("sequence").to_pylist()
        for seq in sequences:
            if not (isinstance(seq, str) and len(seq) > 0):
                raise AssertionError(f"Expected non-empty string, got {seq!r}")

    def test_charge_values_are_valid(self, feature_table):
        charges = feature_table.column("charge").to_pylist()
        for charge in charges:
            if not 1 <= charge <= 10:
                raise AssertionError(f"Charge out of range: {charge}")

    def test_intensities_are_nonnegative(self, feature_table):
        rows = feature_table.column("intensities").to_pylist()
        for row_intensities in rows:
            if row_intensities is None:
                continue
            for entry in row_intensities:
                if entry["intensity"] < 0:
                    raise AssertionError(f"Negative intensity: {entry['intensity']}")

    def test_run_file_names_are_nonempty(self, feature_table):
        run_names = feature_table.column("run_file_name").to_pylist()
        for name in run_names:
            if not (isinstance(name, str) and len(name) > 0):
                raise AssertionError(f"Expected non-empty run name, got {name!r}")

    def test_schema_validation(self, feature_table):
        """Verify that the output conforms to the FeatureSchema."""
        from qpx.core.data import FeatureSchema

        errors = FeatureSchema.validate(feature_table)
        if errors:
            raise AssertionError(f"Schema validation errors: {errors}")

    def test_tmt_nc_isomer_channels_have_distinct_intensities(self, feature_table):
        """Regression test: TMT127N and TMT127C must not systematically share the
        same intensity value.

        The alphabetical sort + i+1 fallback bug in _build_intensities caused
        TMT127C to read from the TMT127N MaxQuant column (and vice-versa),
        producing duplicate intensities for 96 %+ of features in real datasets.
        """
        rows = feature_table.column("intensities").to_pylist()

        equal_pairs = 0
        total_pairs = 0

        for row_intensities in rows:
            if not row_intensities:
                continue
            by_label = {e["label"]: e["intensity"] for e in row_intensities if e["intensity"] > 0}
            n_val = by_label.get("TMT127N")
            c_val = by_label.get("TMT127C")
            if n_val is not None and c_val is not None:
                total_pairs += 1
                if abs(n_val - c_val) < 1e-4:
                    equal_pairs += 1

        if total_pairs == 0:
            pytest.skip("No features with both TMT127N and TMT127C detected")

        duplicate_rate = equal_pairs / total_pairs
        if duplicate_rate >= 0.10:
            raise AssertionError(
                f"TMT127N == TMT127C in {equal_pairs}/{total_pairs} features "
                f"({100 * duplicate_rate:.1f}%) — N/C isomer channel assignment bug detected"
            )


# ---------------------------------------------------------------------------
# fixed_mod_only filter tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestMaxQuantFixedModOnly:
    """Validate that fixed_mod_only=True filters out variable-modification rows."""

    @pytest.fixture(scope="class")
    def fixed_mod_output(self, tmp_path_factory):
        from qpx.converters.maxquant.converter import MaxQuantConverter

        if not _EVIDENCE.exists():
            pytest.skip(f"Test data not found: {EXAMPLES_DIR}")

        output = tmp_path_factory.mktemp("maxquant_fixed_mod")
        converter = MaxQuantConverter()
        converter.convert(
            output_folder=output,
            evidence_file=_EVIDENCE,
            sdrf_file=_SDRF if _SDRF.exists() else None,
            output_prefix="fixed_mod_test",
            fixed_mod_only=True,
        )
        return output

    @pytest.fixture(scope="class")
    def fixed_feature_table(self, fixed_mod_output):
        path = fixed_mod_output / "fixed_mod_test.feature.parquet"
        if not path.exists():
            pytest.skip("fixed_mod feature.parquet was not produced")
        return pq.read_table(str(path))

    def test_fixed_mod_has_fewer_rows_than_unrestricted(self, feature_table, fixed_feature_table):
        """Filtering to fixed mods only must reduce the feature count."""
        if fixed_feature_table.num_rows >= feature_table.num_rows:
            raise AssertionError(
                f"fixed_mod_only produced {fixed_feature_table.num_rows} rows, "
                f"expected fewer than unrestricted {feature_table.num_rows}"
            )

    def test_fixed_mod_contains_no_variable_modifications(self, fixed_feature_table):
        """All features in fixed-mod output must have Unmodified or Carbamidomethyl (C) only.

        Uses startswith("carbamidomethyl") to match ProForma-parsed names regardless
        of whether the exact string is "Carbamidomethyl", "Carbamidomethyl (C)", etc.
        """
        mods_col = fixed_feature_table.column("modifications").to_pylist()
        for mods in mods_col:
            if mods is None:
                continue
            for mod in mods:
                name = (mod.get("modification_name") or "").strip().lower()
                if name:
                    if not name.startswith("carbamidomethyl"):
                        raise AssertionError(f"Variable modification found in fixed-mod output: {mod['modification_name']!r}")


# ---------------------------------------------------------------------------
# Ontology conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestMaxQuantOntologyConversion:
    """Validate ontology.parquet output from MaxQuant conversion.

    The ontology file is only written when scores are discovered during
    conversion. It may legitimately not exist if no scores were tracked.
    """

    def test_ontology_file_exists_or_no_scores(self, converted_output):
        path = converted_output / f"{_PREFIX}.ontology.parquet"
        if not path.exists():
            pytest.skip("ontology.parquet not written (no scores discovered)")

        table = pq.read_table(str(path))
        if table.num_rows == 0:
            raise AssertionError("ontology.parquet exists but is empty")

    def test_ontology_columns(self, converted_output):
        path = converted_output / f"{_PREFIX}.ontology.parquet"
        if not path.exists():
            pytest.skip("ontology.parquet not written (no scores discovered)")

        table = pq.read_table(str(path))
        column_names = set(table.column_names)
        expected = {"field_name", "view"}
        missing = expected - column_names
        if missing:
            raise AssertionError(f"Missing ontology columns: {missing}")

    def test_schema_validation(self, converted_output):
        from qpx.core.data import OntologySchema

        path = converted_output / f"{_PREFIX}.ontology.parquet"
        if not path.exists():
            pytest.skip("ontology.parquet not written (no scores discovered)")
        table = pq.read_table(str(path))
        errors = OntologySchema.validate(table)
        if errors:
            raise AssertionError(f"Schema validation errors: {errors}")
