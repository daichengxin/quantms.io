"""CDAP converter integration tests.

Runs the full CDAP conversion pipeline once (module scope) and validates
that psm.parquet, feature.parquet, and pg.parquet are written with correct
schemas and plausible values.

By default uses the small bundled fixture at ``tests/examples/cdap/``.
Set environment variable ``CDAP_TEST_DATA_DIR`` to a full study directory
(e.g. ``/data/shenyufei/Bigbio_data/CPTAC/stage2/PDC000227``) for
full-scale validation.
"""

import os
from pathlib import Path

import pyarrow.parquet as pq
import pytest

_FIXTURE_DIR = Path(__file__).resolve().parent.parent / "examples" / "cdap"
_PSM_DIR = Path(os.environ.get("CDAP_TEST_DATA_DIR", str(_FIXTURE_DIR)))
_PREFIX = "cdap_test"


# ---------------------------------------------------------------------------
# Module-scoped fixture: run conversion once, share outputs across tests
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def converted_output(tmp_path_factory):
    """Run CDAP conversion once for all tests in this module."""
    from qpx.converters.cdap import CdapConverter

    if not _PSM_DIR.exists() or not list(_PSM_DIR.glob("*.psm")):
        pytest.skip(f"Test data not found: {_PSM_DIR}")

    output = tmp_path_factory.mktemp("cdap_output")

    converter = CdapConverter(max_memory="8GB", max_cpus=24)
    converter.convert(
        psm_dir=_PSM_DIR,
        output_folder=output,
        output_prefix=_PREFIX,
        project_accession="PDC000227",
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
# PSM conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestCdapPsmConversion:
    """Validate psm.parquet output from CDAP conversion."""

    def test_file_exists(self, converted_output):
        assert (converted_output / f"{_PREFIX}.psm.parquet").exists()

    def test_has_rows(self, psm_table):
        assert psm_table.num_rows > 0

    def test_sequence_values_nonempty(self, psm_table):
        for seq in psm_table.column("sequence").to_pylist():
            assert isinstance(seq, str) and len(seq) > 0

    def test_charge_in_range(self, psm_table):
        for charge in psm_table.column("charge").to_pylist():
            assert 1 <= charge <= 15, f"Charge out of range: {charge}"

    def test_calculated_mz_positive(self, psm_table):
        for mz in psm_table.column("calculated_mz").to_pylist():
            if mz is not None:
                assert mz > 0, f"Invalid calculated_mz: {mz}"

    def test_run_names_nonempty(self, psm_table):
        for name in psm_table.column("run_file_name").to_pylist():
            assert isinstance(name, str) and len(name) > 0

    def test_decoy_ratio_reasonable(self, psm_table):
        decoys = sum(1 for x in psm_table.column("is_decoy").to_pylist() if x)
        ratio = decoys / psm_table.num_rows
        assert ratio <= 0.15, f"Decoy ratio too high: {ratio:.2%}"

    def test_schema_validation(self, psm_table):
        from qpx.core.data import PsmSchema

        errors = PsmSchema.validate(psm_table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# Feature conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestCdapFeatureConversion:
    """Validate feature.parquet output from CDAP conversion."""

    def test_file_exists(self, converted_output):
        assert (converted_output / f"{_PREFIX}.feature.parquet").exists()

    def test_has_rows(self, feature_table):
        assert feature_table.num_rows > 0

    def test_sequence_values_nonempty(self, feature_table):
        for seq in feature_table.column("sequence").to_pylist():
            assert isinstance(seq, str) and len(seq) > 0

    def test_intensities_nonnegative(self, feature_table):
        for row_ints in feature_table.column("intensities").to_pylist():
            if row_ints is None:
                continue
            for entry in row_ints:
                assert entry["intensity"] >= 0, f"Negative intensity: {entry}"

    def test_intensities_have_tmt_labels(self, feature_table):
        """At least some features should have TMT10 labels."""
        labels_seen = set()
        for row_ints in feature_table.column("intensities").to_pylist()[:1000]:
            if row_ints is None:
                continue
            for entry in row_ints:
                labels_seen.add(entry["label"])
        assert any("TMT" in lab for lab in labels_seen), f"No TMT labels found: {labels_seen}"

    def test_anchor_protein_nonempty(self, feature_table):
        for ap in feature_table.column("anchor_protein").to_pylist():
            assert isinstance(ap, str) and len(ap) > 0

    def test_schema_validation(self, feature_table):
        from qpx.core.data import FeatureSchema

        errors = FeatureSchema.validate(feature_table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# PG conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestCdapPgConversion:
    """Validate pg.parquet output from CDAP conversion."""

    def test_file_exists(self, converted_output):
        assert (converted_output / f"{_PREFIX}.pg.parquet").exists()

    def test_has_rows(self, pg_table):
        assert pg_table.num_rows > 0

    def test_anchor_protein_nonempty(self, pg_table):
        for ap in pg_table.column("anchor_protein").to_pylist():
            assert isinstance(ap, str) and len(ap) > 0

    def test_intensities_nonnegative(self, pg_table):
        # pg is flattened since 1.1: scalar intensity column (one row per label).
        for intensity in pg_table.column("intensity").to_pylist():
            if intensity is None:
                continue
            assert intensity >= 0, f"Negative intensity: {intensity}"

    def test_pg_accessions_contain_anchor(self, pg_table):
        """The anchor protein must appear in pg_accessions."""
        anchors = pg_table.column("anchor_protein").to_pylist()
        pg_accs = pg_table.column("pg_accessions").to_pylist()
        for anchor, accs in zip(anchors, pg_accs):
            acc_names = [a["accession"] if isinstance(a, dict) else a for a in (accs or [])]
            assert anchor in acc_names, f"Anchor {anchor} not in pg_accessions {acc_names[:5]}"

    def test_schema_validation(self, pg_table):
        from qpx.core.data import PgSchema

        errors = PgSchema.validate(pg_table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# Ontology & metadata tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestCdapMetadata:
    """Validate ontology, provenance, and dataset parquet files."""

    def test_ontology_exists_and_nonempty(self, converted_output):
        path = converted_output / f"{_PREFIX}.ontology.parquet"
        if not path.exists():
            pytest.skip("ontology.parquet not written")
        table = pq.read_table(str(path))
        assert table.num_rows > 0

    def test_provenance_exists(self, converted_output):
        path = converted_output / f"{_PREFIX}.provenance.parquet"
        assert path.exists()
        table = pq.read_table(str(path))
        assert table.num_rows > 0

    def test_provenance_does_not_claim_sample_run(self, converted_output):
        """CDAP does not write sample/run, so provenance must not claim them."""
        path = converted_output / f"{_PREFIX}.provenance.parquet"
        table = pq.read_table(str(path))
        all_views = {view for views in table.column("output_views").to_pylist() for view in (views or [])}
        assert "sample" not in all_views
        assert "run" not in all_views

    def test_dataset_exists(self, converted_output):
        path = converted_output / f"{_PREFIX}.dataset.parquet"
        assert path.exists()
        table = pq.read_table(str(path))
        assert table.num_rows > 0


# ---------------------------------------------------------------------------
# LFQ compatibility tests
# ---------------------------------------------------------------------------


def test_cdap_lfq_precursor_area_uses_lfq_label(tmp_path):
    """CDAP LFQ PrecursorArea must be emitted with an LFQ label."""
    psm_dir = tmp_path / "lfq_psm"
    psm_dir.mkdir()
    psm_path = psm_dir / "sample_lfq.psm"
    psm_path.write_text(
        "\t".join(
            [
                "FileName",
                "ScanNum",
                "QueryPrecursorMz",
                "OriginalPrecursorMz",
                "PrecursorError(ppm)",
                "QueryCharge",
                "OriginalCharge",
                "PrecursorScanNum",
                "PrecursorArea",
                "PrecursorRelAb",
                "RTAtPrecursorHalfElution",
                "PeptideSequence",
                "AmbiguousMatch",
                "Protein",
                "DeNovoScore",
                "MSGFScore",
                "Evalue",
                "Qvalue",
                "PepQvalue",
                "PrecursorPurity",
                "FractionDecomposition",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "run_lfq.raw",
                "1001",
                "500.25",
                "500.25",
                "1.2",
                "2",
                "2",
                "1000?",
                "12345.6",
                "0.1",
                "321.0",
                "PEPTIDEK",
                "0",
                "NP_000001.1(pre=K,post=A)",
                "50",
                "40",
                "0.001",
                "0.0",
                "0.0",
                "99.0,99.0",
                "99.0,99.0",
            ]
        )
        + "\n"
    )

    from qpx.converters.cdap.feature_adapter import CdapFeatureAdapter

    output_path = tmp_path / "lfq.feature.parquet"
    with CdapFeatureAdapter() as adapter:
        adapter.convert(psm_dir=str(psm_dir), output_path=str(output_path))

    table = pq.read_table(str(output_path))
    labels = {entry["label"] for row in table.column("intensities").to_pylist() for entry in (row or [])}
    assert labels == {"LFQ"}


# ---------------------------------------------------------------------------
# Peptidoform unit tests
# ---------------------------------------------------------------------------


class TestCdapPeptidoform:
    """Unit tests for the ProForma conversion logic."""

    def test_nterm_label_strip(self):
        from qpx.converters.cdap.peptidoform import strip_label_tags

        stripped = strip_label_tags("+229.163PEPTIDEK+229.163")
        assert stripped == "PEPTIDEK"

    def test_internal_ptm(self):
        from qpx.converters.cdap.peptidoform import cdap_to_proforma

        result = cdap_to_proforma("+229.163PEPTM+15.995IDEK+229.163")
        assert "UNIMOD:" in result or "+15.995" in result

    def test_no_labels(self):
        from qpx.converters.cdap.peptidoform import strip_label_tags

        stripped = strip_label_tags("PEPTIDEK")
        assert stripped == "PEPTIDEK"

    def test_plain_sequence_roundtrip(self):
        from qpx.converters.cdap.peptidoform import cdap_to_proforma

        result = cdap_to_proforma("PEPTIDEK")
        assert result == "PEPTIDEK"
