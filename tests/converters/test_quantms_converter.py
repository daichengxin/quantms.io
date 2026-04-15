"""QuantMS converter integration tests using DDA-LFQ full dataset.

Two tiers of tests:

- **SDRF-only tests** (sample, run, ontology) run quickly (~2s) and are always
  included in the regular test suite.
- **mzTab-based tests** (psm, feature, pg) require parsing a 41 MB gzipped
  mzTab file and are marked ``large_data``.  Run them explicitly::

      python -m pytest tests/converters/test_quantms_converter.py -v -m "large_data"
"""

from pathlib import Path

import pyarrow.parquet as pq
import pytest

EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "examples" / "quantms" / "dda-lfq-full"

_MZTAB = EXAMPLES_DIR / "PXD007683-LFQ.sdrf_openms_design_openms.mzTab.gz"
_SDRF = EXAMPLES_DIR / "PXD007683-LFQ.sdrf.tsv"
_MSSTATS = EXAMPLES_DIR / "PXD007683-LFQ.sdrf_openms_design_msstats_in.csv.gz"

_PREFIX = "quantms_test"


# ---------------------------------------------------------------------------
# Fast fixture: SDRF conversion only (sample + run + ontology)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def sdrf_output(tmp_path_factory):
    """Convert SDRF only — fast (~2s), produces sample + run + ontology."""
    from qpx.converters.quantms.converter import QuantMSConverter

    if not _SDRF.exists():
        pytest.skip(f"Test data not found: {_SDRF}")
    if not _MZTAB.exists():
        pytest.skip(f"Test data not found: {_MZTAB}")

    output = tmp_path_factory.mktemp("quantms_sdrf_output")

    converter = QuantMSConverter(
        mztab_path=str(_MZTAB),
        sdrf_file=str(_SDRF),
    )
    converter.convert(
        output_folder=output,
        output_prefix=_PREFIX,
        structures=[],  # SDRF is always converted
    )

    return output


@pytest.fixture(scope="module")
def sample_table(sdrf_output):
    path = sdrf_output / f"{_PREFIX}.sample.parquet"
    if not path.exists():
        pytest.skip("sample.parquet was not produced")
    return pq.read_table(str(path))


@pytest.fixture(scope="module")
def run_table(sdrf_output):
    path = sdrf_output / f"{_PREFIX}.run.parquet"
    if not path.exists():
        pytest.skip("run.parquet was not produced")
    return pq.read_table(str(path))


# ---------------------------------------------------------------------------
# Slow fixture: full mzTab conversion (psm + feature + pg)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def converted_output(tmp_path_factory):
    """Run full QuantMS conversion — slow (~27 min for 41 MB mzTab).

    Only triggered by tests marked ``large_data``.
    """
    from qpx.converters.quantms.converter import QuantMSConverter

    if not _MZTAB.exists():
        pytest.skip(f"Test data not found: {_MZTAB}")
    if not _SDRF.exists():
        pytest.skip(f"Test data not found: {_SDRF}")

    output = tmp_path_factory.mktemp("quantms_full_output")

    converter = QuantMSConverter(
        mztab_path=str(_MZTAB),
        sdrf_file=str(_SDRF),
        msstats_file=str(_MSSTATS) if _MSSTATS.exists() else None,
    )
    converter.convert(
        output_folder=output,
        output_prefix=_PREFIX,
        structures=["psm", "feature", "pg"],
    )

    return output


@pytest.fixture(scope="module")
def psm_table(converted_output):
    path = converted_output / f"{_PREFIX}.psm.parquet"
    if not path.exists():
        pytest.skip("psm.parquet was not produced")
    return pq.read_table(str(path))


@pytest.fixture(scope="module")
def feature_table(converted_output):
    path = converted_output / f"{_PREFIX}.feature.parquet"
    if not path.exists():
        pytest.skip("feature.parquet was not produced")
    return pq.read_table(str(path))


@pytest.fixture(scope="module")
def pg_table(converted_output):
    path = converted_output / f"{_PREFIX}.pg.parquet"
    if not path.exists():
        pytest.skip("pg.parquet was not produced")
    return pq.read_table(str(path))


# ---------------------------------------------------------------------------
# Sample conversion tests (fast — uses sdrf_output fixture)
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestQuantMSSampleConversion:
    """Validate sample.parquet output from QuantMS SDRF conversion."""

    def test_file_exists(self, sdrf_output):
        path = sdrf_output / f"{_PREFIX}.sample.parquet"
        assert path.exists(), "sample.parquet was not created"

    def test_has_rows(self, sample_table):
        assert sample_table.num_rows > 0, "sample.parquet is empty"

    def test_schema_validation(self, sample_table):
        from qpx.core.data import SampleSchema

        errors = SampleSchema.validate(sample_table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# Run conversion tests (fast — uses sdrf_output fixture)
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestQuantMSRunConversion:
    """Validate run.parquet output from QuantMS SDRF conversion."""

    def test_file_exists(self, sdrf_output):
        path = sdrf_output / f"{_PREFIX}.run.parquet"
        assert path.exists(), "run.parquet was not created"

    def test_has_rows(self, run_table):
        assert run_table.num_rows > 0, "run.parquet is empty"

    def test_schema_validation(self, run_table):
        from qpx.core.data import RunSchema

        errors = RunSchema.validate(run_table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# Ontology conversion tests (fast — uses sdrf_output fixture)
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestQuantMSOntologyConversion:
    """Validate ontology.parquet from SDRF conversion."""

    def test_ontology_file_exists(self, sdrf_output):
        path = sdrf_output / f"{_PREFIX}.ontology.parquet"
        if not path.exists():
            pytest.skip("ontology.parquet not written (no scores discovered)")
        table = pq.read_table(str(path))
        assert table.num_rows > 0, "ontology.parquet exists but is empty"

    def test_ontology_columns(self, sdrf_output):
        path = sdrf_output / f"{_PREFIX}.ontology.parquet"
        if not path.exists():
            pytest.skip("ontology.parquet not written (no scores discovered)")
        table = pq.read_table(str(path))
        column_names = set(table.column_names)
        expected = {"field_name", "view"}
        missing = expected - column_names
        assert not missing, f"Missing ontology columns: {missing}"

    def test_schema_validation(self, sdrf_output):
        from qpx.core.data import OntologySchema

        path = sdrf_output / f"{_PREFIX}.ontology.parquet"
        if not path.exists():
            pytest.skip("ontology.parquet not written (no scores discovered)")
        table = pq.read_table(str(path))
        errors = OntologySchema.validate(table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# PSM conversion tests (slow — uses converted_output fixture)
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestQuantMSPsmConversion:
    """Validate psm.parquet output from QuantMS conversion."""

    def test_file_exists(self, converted_output):
        path = converted_output / f"{_PREFIX}.psm.parquet"
        assert path.exists(), "psm.parquet was not created"

    def test_has_rows(self, psm_table):
        assert psm_table.num_rows > 0, "psm.parquet is empty"

    def test_key_columns_present(self, psm_table):
        column_names = set(psm_table.column_names)
        expected = {"sequence", "charge", "calculated_mz", "run_file_name"}
        missing = expected - column_names
        assert not missing, f"Missing columns: {missing}"

    def test_sequence_values_are_nonempty(self, psm_table):
        sequences = psm_table.column("sequence").to_pylist()
        for seq in sequences:
            assert isinstance(seq, str) and len(seq) > 0, f"Expected non-empty string, got {seq!r}"

    def test_charge_values_are_valid(self, psm_table):
        charges = psm_table.column("charge").to_pylist()
        for charge in charges:
            assert 1 <= charge <= 10, f"Charge out of range: {charge}"

    def test_calculated_mz_values_are_nonnegative(self, psm_table):
        mz_values = psm_table.column("calculated_mz").to_pylist()
        for mz in mz_values:
            assert mz is not None and mz >= 0, f"Invalid calculated_mz value: {mz}"

    def test_run_file_names_are_nonempty(self, psm_table):
        run_names = psm_table.column("run_file_name").to_pylist()
        for name in run_names:
            assert isinstance(name, str) and len(name) > 0

    def test_schema_validation(self, psm_table):
        from qpx.core.data import PsmSchema

        errors = PsmSchema.validate(psm_table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# Feature conversion tests (slow — uses converted_output fixture)
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestQuantMSFeatureConversion:
    """Validate feature.parquet output from QuantMS conversion."""

    def test_file_exists(self, converted_output):
        path = converted_output / f"{_PREFIX}.feature.parquet"
        assert path.exists(), "feature.parquet was not created"

    def test_has_rows(self, feature_table):
        assert feature_table.num_rows > 0, "feature.parquet is empty"

    def test_key_columns_present(self, feature_table):
        column_names = set(feature_table.column_names)
        expected = {"sequence", "charge", "calculated_mz", "run_file_name"}
        missing = expected - column_names
        assert not missing, f"Missing columns: {missing}"

    def test_sequence_values_are_nonempty(self, feature_table):
        sequences = feature_table.column("sequence").to_pylist()
        for seq in sequences:
            assert isinstance(seq, str) and len(seq) > 0, f"Expected non-empty string, got {seq!r}"

    def test_run_file_names_are_nonempty(self, feature_table):
        run_names = feature_table.column("run_file_name").to_pylist()
        for name in run_names:
            assert isinstance(name, str) and len(name) > 0

    def test_schema_validation(self, feature_table):
        from qpx.core.data import FeatureSchema

        errors = FeatureSchema.validate(feature_table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# Protein group conversion tests (slow — uses converted_output fixture)
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestQuantMSPgConversion:
    """Validate pg.parquet output from QuantMS conversion."""

    def test_file_exists(self, converted_output):
        path = converted_output / f"{_PREFIX}.pg.parquet"
        assert path.exists(), "pg.parquet was not created"

    def test_has_rows(self, pg_table):
        assert pg_table.num_rows > 0, "pg.parquet is empty"

    def test_key_columns_present(self, pg_table):
        column_names = set(pg_table.column_names)
        expected = {"anchor_protein", "peptides"}
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

    def test_schema_validation(self, pg_table):
        from qpx.core.data import PgSchema

        errors = PgSchema.validate(pg_table)
        assert not errors, f"Schema validation errors: {errors}"
