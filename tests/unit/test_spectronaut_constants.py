"""Unit tests for Spectronaut PTM constants."""

from qpx.converters.spectronaut.constants import to_modifications, to_proforma


class TestToProforma:
    """Test to_proforma conversion from Spectronaut format."""

    def test_unmodified(self):
        assert to_proforma("_EVVEAHVDQK_") == "EVVEAHVDQK"

    def test_empty(self):
        assert to_proforma("") == ""

    def test_single_carbamidomethyl(self):
        result = to_proforma("_YQC[Carbamidomethyl (C)]VVLTEMK_")
        assert result == "YQC[UNIMOD:4]VVLTEMK"

    def test_double_carbamidomethyl(self):
        result = to_proforma("_PIGLC[Carbamidomethyl (C)]C[Carbamidomethyl (C)]IAPVLAAK_")
        assert result == "PIGLC[UNIMOD:4]C[UNIMOD:4]IAPVLAAK"

    def test_oxidation(self):
        result = to_proforma("_GM[Oxidation (M)]ITVTDPDLIEK_")
        assert result == "GM[UNIMOD:35]ITVTDPDLIEK"

    def test_nterm_acetyl(self):
        result = to_proforma("_[Acetyl (Protein N-term)]PEPTIDEK_")
        assert result == "[UNIMOD:1]-PEPTIDEK"

    def test_nterm_plus_internal_mod(self):
        result = to_proforma("_[Acetyl (Protein N-term)]M[Oxidation (M)]PEPTIDEK_")
        assert result == "[UNIMOD:1]-M[UNIMOD:35]PEPTIDEK"

    def test_unknown_mod_passes_through(self):
        result = to_proforma("_PEP[SomeUnknownMod]TIDEK_")
        assert "SomeUnknownMod" in result

    def test_caching(self):
        a = to_proforma("_EVVEAHVDQK_")
        b = to_proforma("_EVVEAHVDQK_")
        assert a == b


class TestToModifications:
    """Test to_modifications returning QPX modification dicts."""

    def test_unmodified_returns_none(self):
        pf, mods = to_modifications("_PEPTIDEK_", "PEPTIDEK")
        assert pf == "PEPTIDEK"
        assert mods is None

    def test_carbamidomethyl_returns_mod(self):
        pf, mods = to_modifications("_YQC[Carbamidomethyl (C)]VVLTEMK_", "YQCVVLTEMK")
        assert "UNIMOD:4" in pf
        assert mods is not None
        assert len(mods) >= 1
        assert any("UNIMOD:4" in str(m.get("accession", "")) for m in mods)

    def test_nterm_acetyl(self):
        pf, mods = to_modifications("_[Acetyl (Protein N-term)]PEPTIDEK_", "PEPTIDEK")
        assert "[UNIMOD:1]-" in pf
        assert mods is not None
        positions = []
        for mod in mods:
            for pos in mod.get("positions", []):
                positions.append(pos["position"])
        assert 0 in positions
