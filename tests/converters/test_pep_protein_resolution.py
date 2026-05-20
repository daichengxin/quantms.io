"""Unit tests for PEP-section peptide→protein resolution logic.

Tests cover:
- ``_build_peptide_protein_map()`` with a real peptides table
- ``_build_peptide_protein_map()`` with a fallback (empty / missing columns)
- protein-group resolution in ``_rows_to_feature_records()``
"""

import duckdb
import pytest

from qpx.converters.quantms.feature_adapter import QuantmsFeatureAdapter


@pytest.fixture()
def adapter():
    """Create an adapter with an in-memory DuckDB connection."""
    conn = duckdb.connect(":memory:")
    a = QuantmsFeatureAdapter(conn=conn)
    return a


class TestBuildPeptideProteinMap:
    """Tests for _build_peptide_protein_map()."""

    def test_unambiguous_mapping(self, adapter):
        """Peptides mapped to a single protein return a correct dict."""
        adapter._conn.execute("""
            CREATE TABLE peptides (sequence VARCHAR, accession VARCHAR)
        """)
        adapter._conn.execute("""
            INSERT INTO peptides VALUES
                ('PEPTIDEK', 'P12345'),
                ('ANOTHERK', 'P67890'),
                ('PEPTIDEK', 'P12345')
        """)
        result = adapter._build_peptide_protein_map()
        assert result == {"PEPTIDEK": "P12345", "ANOTHERK": "P67890"}

    def test_shared_peptide_excluded(self, adapter):
        """Peptides mapping to multiple proteins are excluded."""
        adapter._conn.execute("""
            CREATE TABLE peptides (sequence VARCHAR, accession VARCHAR)
        """)
        adapter._conn.execute("""
            INSERT INTO peptides VALUES
                ('SHAREDK', 'P11111'),
                ('SHAREDK', 'P22222'),
                ('UNIQUEK', 'P33333')
        """)
        result = adapter._build_peptide_protein_map()
        assert "SHAREDK" not in result
        assert result == {"UNIQUEK": "P33333"}

    def test_no_peptides_table(self, adapter):
        """Returns empty dict when peptides table does not exist."""
        result = adapter._build_peptide_protein_map()
        assert result == {}

    def test_fallback_table_missing_columns(self, adapter):
        """Returns empty dict when peptides table lacks sequence/accession."""
        adapter._conn.execute("CREATE TABLE peptides (id INTEGER)")
        result = adapter._build_peptide_protein_map()
        assert result == {}


class TestProteinGroupResolution:
    """Tests for protein-group resolution in _rows_to_feature_records()."""

    def _make_row(self, protein_name="P11;P22", sequence="TESTPEPTIDE"):
        """Build a minimal row tuple matching the SQL query column order.

        Column order:
        0: peptidoform, 1: sequence, 2: run_file_name, 3: charge,
        4: intensity, 5: msstats_rt, 6: protein_name,
        7-13: PSM fields, 14: pg_global_qvalue, 15: gene_name,
        16: modifications_json
        """
        return (
            sequence,  # 0: peptidoform
            sequence,  # 1: sequence
            "run1",  # 2: run_file_name
            2,  # 3: charge
            1000.0,  # 4: intensity
            12.5,  # 5: msstats_rt
            protein_name,  # 6: protein_name
            None,  # 7: psm_pep
            500.0,  # 8: psm_calc_mz
            500.1,  # 9: psm_obs_mz
            False,  # 10: psm_is_decoy
            None,  # 11: psm_scan
            None,  # 12: psm_id_run
            None,  # 13: psm_rt
            0.01,  # 14: pg_global_qvalue
            "GN1",  # 15: gene_name
            None,  # 16: modifications_json
        )

    def test_resolves_ambiguous_group(self, adapter):
        """Multi-accession protein_name 'A;B' is resolved to 'A' via PEP map."""
        adapter._pep_protein_map = {"TESTPEPTIDE": "P11"}
        adapter._enzyme_name = None

        rows = [self._make_row(protein_name="P11;P22", sequence="TESTPEPTIDE")]
        records = adapter._rows_to_feature_records(rows)

        assert len(records) == 1
        accessions = [a["accession"] for a in records[0]["pg_accessions"]]
        assert accessions == ["P11"]
        assert records[0]["anchor_protein"] == "P11"
        assert records[0]["unique"] is True

    def test_keeps_single_accession(self, adapter):
        """Single-accession protein_name is unchanged."""
        adapter._pep_protein_map = {"TESTPEPTIDE": "P11"}
        adapter._enzyme_name = None

        rows = [self._make_row(protein_name="P11", sequence="TESTPEPTIDE")]
        records = adapter._rows_to_feature_records(rows)

        assert len(records) == 1
        accessions = [a["accession"] for a in records[0]["pg_accessions"]]
        assert accessions == ["P11"]

    def test_no_map_keeps_group(self, adapter):
        """Without PEP map, multi-accession group is preserved."""
        adapter._pep_protein_map = {}
        adapter._enzyme_name = None

        rows = [self._make_row(protein_name="P11;P22", sequence="TESTPEPTIDE")]
        records = adapter._rows_to_feature_records(rows)

        assert len(records) == 1
        accessions = [a["accession"] for a in records[0]["pg_accessions"]]
        assert accessions == ["P11", "P22"]
        assert records[0]["unique"] is False
