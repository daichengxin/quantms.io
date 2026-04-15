"""Unit tests for converter utility methods added in recent commits.

Tests cover:
- BaseConverter._escape_path (SQL path injection prevention)
- MaxQuantPgAdapter._extract_protein_names / _extract_gene_names
- MaxQuantFeatureAdapter._detect_pg_columns / _extract_qvalue_map / _extract_gene_map
"""

import pandas as pd

from qpx.converters.base import BaseConverter
from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter
from qpx.converters.maxquant.pg_adapter import MaxQuantPgAdapter

# ---------------------------------------------------------------------------
# BaseConverter._escape_path
# ---------------------------------------------------------------------------


class TestEscapePath:
    """Validate SQL path escaping."""

    def test_normal_path(self):
        result = BaseConverter._escape_path("/data/results.txt")
        if result != "/data/results.txt":
            raise AssertionError(f"Expected '/data/results.txt', got {result!r}")

    def test_path_with_single_quote(self):
        result = BaseConverter._escape_path("O'Brien/data.txt")
        if result != "O''Brien/data.txt":
            raise AssertionError(f"Expected \"O''Brien/data.txt\", got {result!r}")

    def test_path_with_multiple_quotes(self):
        result = BaseConverter._escape_path("a'b'c.txt")
        if result != "a''b''c.txt":
            raise AssertionError(f"Expected \"a''b''c.txt\", got {result!r}")

    def test_windows_path(self):
        result = BaseConverter._escape_path("C:\\Users\\test\\file.txt")
        if result != "C:\\Users\\test\\file.txt":
            raise AssertionError(f"Expected windows path unchanged, got {result!r}")

    def test_empty_path(self):
        result = BaseConverter._escape_path("")
        if result != "":
            raise AssertionError(f"Expected empty string, got {result!r}")


# ---------------------------------------------------------------------------
# MaxQuantPgAdapter._extract_gene_names
# ---------------------------------------------------------------------------


class TestExtractGeneNames:
    """Validate GN= tag extraction from FASTA headers."""

    def test_single_gene(self):
        result = MaxQuantPgAdapter._extract_gene_names("sp|P12345|PROT_HUMAN Protein X OS=Homo sapiens GN=TP53 PE=1 SV=1")
        if result != ["TP53"]:
            raise AssertionError(f"Unexpected result: {result!r}")

    def test_multiple_headers(self):
        result = MaxQuantPgAdapter._extract_gene_names("sp|P12345|X GN=TP53 PE=1;sp|Q99999|Y GN=BRCA1 PE=1")
        if result != ["TP53", "BRCA1"]:
            raise AssertionError(f"Unexpected result: {result!r}")

    def test_no_gene_name(self):
        result = MaxQuantPgAdapter._extract_gene_names("sp|P12345|PROT_HUMAN Protein X OS=Homo sapiens PE=1 SV=1")
        if result is not None:
            raise AssertionError(f"Expected None, got {result!r}")

    def test_duplicate_gene_names(self):
        result = MaxQuantPgAdapter._extract_gene_names("sp|P12345|X GN=TP53 PE=1;sp|Q99999|Y GN=TP53 PE=1")
        if result != ["TP53"]:
            raise AssertionError(f"Unexpected result: {result!r}")

    def test_empty_string(self):
        result = MaxQuantPgAdapter._extract_gene_names("")
        if result is not None:
            raise AssertionError(f"Expected None, got {result!r}")


# ---------------------------------------------------------------------------
# MaxQuantFeatureAdapter._detect_pg_columns
# ---------------------------------------------------------------------------


class TestDetectPgColumns:
    """Validate column detection in proteinGroups DataFrame."""

    def test_standard_columns(self):
        df = pd.DataFrame({"Protein IDs": ["P1"], "Q-value": [0.01], "Gene names": ["TP53"]})
        acc, qval, gene = MaxQuantFeatureAdapter._detect_pg_columns(df)
        if acc != "Protein IDs":
            raise AssertionError(f"Expected 'Protein IDs', got {acc!r}")
        if qval != "Q-value":
            raise AssertionError(f"Expected 'Q-value', got {qval!r}")
        if gene != "Gene names":
            raise AssertionError(f"Expected 'Gene names', got {gene!r}")

    def test_alternative_columns(self):
        df = pd.DataFrame({"Majority protein IDs": ["P1"], "q-value": [0.01], "Gene Names": ["TP53"]})
        acc, qval, gene = MaxQuantFeatureAdapter._detect_pg_columns(df)
        if acc != "Majority protein IDs":
            raise AssertionError(f"Expected 'Majority protein IDs', got {acc!r}")
        if qval != "q-value":
            raise AssertionError(f"Expected 'q-value', got {qval!r}")
        if gene != "Gene Names":
            raise AssertionError(f"Expected 'Gene Names', got {gene!r}")

    def test_missing_columns(self):
        df = pd.DataFrame({"some_col": [1]})
        acc, qval, gene = MaxQuantFeatureAdapter._detect_pg_columns(df)
        if acc is not None:
            raise AssertionError(f"Expected None for acc, got {acc!r}")
        if qval is not None:
            raise AssertionError(f"Expected None for qval, got {qval!r}")
        if gene is not None:
            raise AssertionError(f"Expected None for gene, got {gene!r}")


# ---------------------------------------------------------------------------
# MaxQuantFeatureAdapter._extract_qvalue_map
# ---------------------------------------------------------------------------


class TestExtractQvalueMap:
    """Validate q-value map extraction."""

    def test_basic_qvalue_map(self):
        df = pd.DataFrame({"Protein IDs": ["P1;P2", "P3"], "Q-value": [0.01, 0.05]})
        result = MaxQuantFeatureAdapter._extract_qvalue_map(df, "Protein IDs", "Q-value")
        if result["P1"] != 0.01:
            raise AssertionError(f"Expected 0.01 for P1, got {result['P1']!r}")
        if result["P2"] != 0.01:
            raise AssertionError(f"Expected 0.01 for P2, got {result['P2']!r}")
        if result["P3"] != 0.05:
            raise AssertionError(f"Expected 0.05 for P3, got {result['P3']!r}")

    def test_first_occurrence_wins(self):
        df = pd.DataFrame({"Protein IDs": ["P1", "P1"], "Q-value": [0.01, 0.99]})
        result = MaxQuantFeatureAdapter._extract_qvalue_map(df, "Protein IDs", "Q-value")
        if result["P1"] != 0.01:
            raise AssertionError(f"Expected 0.01 for P1, got {result['P1']!r}")

    def test_no_qval_col(self):
        df = pd.DataFrame({"Protein IDs": ["P1"]})
        result = MaxQuantFeatureAdapter._extract_qvalue_map(df, "Protein IDs", None)
        if result != {}:
            raise AssertionError(f"Expected empty dict, got {result!r}")


# ---------------------------------------------------------------------------
# MaxQuantFeatureAdapter._extract_gene_map
# ---------------------------------------------------------------------------


class TestExtractGeneMap:
    """Validate gene map extraction with Fasta header fallback."""

    def test_from_gene_column(self):
        df = pd.DataFrame({"Protein IDs": ["P1;P2"], "Gene names": ["TP53;BRCA1"]})
        result = MaxQuantFeatureAdapter._extract_gene_map(df, "Protein IDs", "Gene names")
        if result["P1"] != ["TP53", "BRCA1"]:
            raise AssertionError(f"Unexpected P1: {result['P1']!r}")
        if result["P2"] != ["TP53", "BRCA1"]:
            raise AssertionError(f"Unexpected P2: {result['P2']!r}")

    def test_fasta_fallback(self):
        df = pd.DataFrame(
            {
                "Protein IDs": ["P1"],
                "Gene names": [None],
                "Fasta headers": ["sp|P1|X GN=TP53 PE=1"],
            }
        )
        result = MaxQuantFeatureAdapter._extract_gene_map(df, "Protein IDs", "Gene names")
        if result["P1"] != ["TP53"]:
            raise AssertionError(f"Unexpected P1: {result['P1']!r}")

    def test_no_gene_info(self):
        df = pd.DataFrame({"Protein IDs": ["P1"], "Gene names": [None]})
        result = MaxQuantFeatureAdapter._extract_gene_map(df, "Protein IDs", "Gene names")
        if result != {}:
            raise AssertionError(f"Expected empty dict, got {result!r}")
