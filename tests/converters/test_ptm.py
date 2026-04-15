"""Test shared PTM utilities (flat ptm.py)."""

from qpx.converters.ptm import (
    UNIMOD_MASS,
    _normalize_peptidoform,
    build_proforma,
    compute_precursor_mz,
    from_proforma,
    mass_to_unimod,
)


class TestUnimodMassRegistry:
    def test_oxidation_mass(self):
        assert abs(UNIMOD_MASS[35] - 15.994915) < 0.0001

    def test_acetyl_mass(self):
        assert abs(UNIMOD_MASS[1] - 42.010565) < 0.0001

    def test_mass_to_unimod_oxidation(self):
        assert mass_to_unimod(15.9949) == 35

    def test_mass_to_unimod_unknown(self):
        assert mass_to_unimod(999.0) is None


class TestBuildProforma:
    def test_no_mods(self):
        assert build_proforma("PEPTIDEK", []) == "PEPTIDEK"

    def test_internal_mod(self):
        assert build_proforma("PEPTMIDEK", [(5, "UNIMOD:35")]) == "PEPTM[UNIMOD:35]IDEK"

    def test_nterm_mod(self):
        assert build_proforma("PEPTIDEK", [(0, "UNIMOD:1")]) == "[UNIMOD:1]-PEPTIDEK"

    def test_empty_sequence(self):
        assert build_proforma("", []) == ""


class TestFromProforma:
    def test_unimod_tag(self):
        meta = {"UNIMOD:35": ("Oxidation", ["M"], ["Anywhere"])}
        result = from_proforma("M[UNIMOD:35]PEPTIDEK", "MPEPTIDEK", meta=meta)
        assert result is not None
        assert len(result) == 1
        assert result[0]["name"] == "Oxidation"
        assert result[0]["accession"] == "UNIMOD:35"

    def test_no_mods(self):
        assert from_proforma("PEPTIDEK", "PEPTIDEK", meta=None) is None

    def test_nterm(self):
        meta = {"UNIMOD:1": ("Acetyl", ["X"], ["N-term"])}
        result = from_proforma("[UNIMOD:1]-PEPTIDEK", "PEPTIDEK", meta=meta)
        assert result is not None
        assert result[0]["positions"][0]["position"] == 0

    def test_site_scores_attached(self):
        """Site scores are attached to the correct modification position."""
        site_scores = {
            1: [
                {
                    "score_name": "phosphors_site_probability",
                    "score_value": 0.95,
                    "higher_better": True,
                }
            ],
        }
        result = from_proforma(
            "M[UNIMOD:35]PEPTIDEK",
            "MPEPTIDEK",
            meta=None,
            site_scores=site_scores,
        )
        assert result is not None
        pos = result[0]["positions"][0]
        assert pos["position"] == 1
        assert pos["scores"] is not None
        assert pos["scores"][0]["score_value"] == 0.95

    def test_site_scores_none_when_not_provided(self):
        """Without site_scores, positions still have scores=None."""
        result = from_proforma("M[UNIMOD:35]PEPTIDEK", "MPEPTIDEK", meta=None)
        assert result is not None
        assert result[0]["positions"][0]["scores"] is None

    def test_site_scores_partial_positions(self):
        """Only positions in the site_scores dict get scores attached."""
        site_scores = {
            5: [
                {
                    "score_name": "phospho_prob",
                    "score_value": 0.8,
                    "higher_better": True,
                }
            ],
        }
        # Two mods: position 1 and position 5; only position 5 has scores
        result = from_proforma(
            "M[UNIMOD:35]PEPTS[UNIMOD:21]IDEK",
            "MPEPTSIDEK",
            meta=None,
            site_scores=site_scores,
        )
        assert result is not None
        # Find the phospho mod (position 5)
        for mod in result:
            for pos in mod["positions"]:
                if pos["position"] == 5:
                    assert pos["scores"] is not None
                    assert pos["scores"][0]["score_value"] == 0.8
                elif pos["position"] == 1:
                    assert pos["scores"] is None


class TestNormalizePeptidoform:
    """Validate mzTab parenthetical -> ProForma bracket conversion."""

    def test_no_parens(self):
        result = _normalize_peptidoform("PEPTIDEK")
        if result != "PEPTIDEK":
            raise AssertionError(f"Unexpected: {result!r}")

    def test_simple_mod(self):
        result = _normalize_peptidoform("C(Carbamidomethyl)R")
        if result != "C[Carbamidomethyl]R":
            raise AssertionError(f"Unexpected: {result!r}")

    def test_nested_parens(self):
        result = _normalize_peptidoform("C(Carbamidomethyl (C))R")
        if result != "C[Carbamidomethyl (C)]R":
            raise AssertionError(f"Unexpected: {result!r}")

    def test_nterm_mod(self):
        result = _normalize_peptidoform("(Acetyl)PEPTIDEK")
        if result != "[Acetyl]PEPTIDEK":
            raise AssertionError(f"Unexpected: {result!r}")

    def test_multiple_mods(self):
        result = _normalize_peptidoform("(Acetyl)M(Oxidation)PEPTC(Carbamidomethyl)K")
        if result != "[Acetyl]M[Oxidation]PEPTC[Carbamidomethyl]K":
            raise AssertionError(f"Unexpected: {result!r}")

    def test_unmatched_open_paren(self):
        result = _normalize_peptidoform("C(Broken")
        if result != "C(Broken":
            raise AssertionError(f"Unexpected: {result!r}")

    def test_empty_string(self):
        result = _normalize_peptidoform("")
        if result != "":
            raise AssertionError(f"Unexpected: {result!r}")


class TestFromProformaMzTab:
    """Validate from_proforma with mzTab parenthetical input."""

    def test_mztab_simple(self):
        meta = {"UNIMOD:4": ("Carbamidomethyl", ["C"], ["Anywhere"])}
        result = from_proforma("C(UNIMOD:4)PEPTIDEK", "CPEPTIDEK", meta=meta)
        if result is None:
            raise AssertionError("Expected mods, got None")
        if result[0]["accession"] != "UNIMOD:4":
            raise AssertionError(f"Unexpected accession: {result[0]['accession']!r}")

    def test_mztab_nterm(self):
        meta = {"UNIMOD:1": ("Acetyl", ["X"], ["N-term"])}
        result = from_proforma("(UNIMOD:1)-PEPTIDEK", "PEPTIDEK", meta=meta)
        if result is None:
            raise AssertionError("Expected mods, got None")
        if result[0]["positions"][0]["position"] != 0:
            raise AssertionError(f"Unexpected position: {result[0]['positions'][0]['position']}")

    def test_mztab_nested_parens(self):
        # Normalization preserves inner parens as-is
        result = from_proforma("C(Carbamidomethyl (C))PEPTIDEK", "CPEPTIDEK", meta=None)
        if result is None:
            raise AssertionError("Expected mods, got None")
        mod_name = result[0]["name"]
        if "Carbamidomethyl (C)" not in mod_name:
            raise AssertionError(f"Inner parens lost: {mod_name!r}")


class TestComputePrecursorMz:
    def test_basic(self):
        mz = compute_precursor_mz("PEPTM(UniMod:35)IDE", 2)
        assert mz is not None
        assert mz > 0

    def test_empty(self):
        assert compute_precursor_mz("", 2) is None
        assert compute_precursor_mz("PEPTIDE", 0) is None
