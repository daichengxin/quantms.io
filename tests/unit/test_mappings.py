"""Tests for the central column-mapping loader."""

import pytest


class TestGetFieldMappings:
    """get_field_mappings returns correct column mapping dicts."""

    def test_diann_feature_returns_dict_of_lists(self):
        from qpx.converters.mappings import get_field_mappings

        result = get_field_mappings("diann", "feature")
        assert isinstance(result, dict)
        for qpx_field, candidates in result.items():
            assert isinstance(candidates, list), f"{qpx_field} must map to a list"
            assert len(candidates) > 0

    def test_diann_feature_has_key_fields(self):
        from qpx.converters.mappings import get_field_mappings

        feature = get_field_mappings("diann", "feature")
        assert "intensity" in feature
        assert feature["intensity"] == ["Precursor.Quantity"]
        assert "rt" in feature
        assert "pg_accessions" in feature

    def test_fragpipe_psm(self):
        from qpx.converters.mappings import get_field_mappings

        psm = get_field_mappings("fragpipe", "psm")
        assert psm["sequence"] == ["Peptide"]
        assert psm["observed_mz"] == ["Observed M/Z"]

    def test_maxquant_pg(self):
        from qpx.converters.mappings import get_field_mappings

        pg = get_field_mappings("maxquant", "pg")
        assert pg["pg_accessions"] == ["Protein IDs"]
        assert "andromeda_score" in pg

    def test_quantms_feature(self):
        from qpx.converters.mappings import get_field_mappings

        feature = get_field_mappings("quantms", "feature")
        assert feature["peptidoform"] == ["PeptideSequence", "peptidoform", "Peptide"]

    def test_unknown_tool_raises(self):
        from qpx.converters.mappings import get_field_mappings

        with pytest.raises(KeyError):
            get_field_mappings("nonexistent", "feature")

    def test_unknown_view_raises(self):
        from qpx.converters.mappings import get_field_mappings

        with pytest.raises(KeyError):
            get_field_mappings("diann", "nonexistent")


class TestGetToolMeta:
    """get_tool_meta returns tool name, versions, and rt_unit."""

    def test_diann_meta(self):
        from qpx.converters.mappings import get_tool_meta

        meta = get_tool_meta("diann")
        assert meta["tool_name"] == "DIA-NN"
        assert meta["tool_versions"] == "1.8+"
        assert meta["rt_unit"] == "minute"

    def test_maxquant_meta(self):
        from qpx.converters.mappings import get_tool_meta

        meta = get_tool_meta("maxquant")
        assert meta["tool_name"] == "MaxQuant"
        assert meta["tool_versions"] == "2.x"

    def test_unknown_tool_raises(self):
        from qpx.converters.mappings import get_tool_meta

        with pytest.raises(KeyError):
            get_tool_meta("nonexistent")


class TestGetExtra:
    """get_extra returns tool-specific config sections."""

    def test_quantms_phospho_site_columns(self):
        from qpx.converters.mappings import get_extra

        phospho = get_extra("quantms", "phospho_site_columns")
        assert isinstance(phospho, dict)
        assert phospho["opt_global_phosphors_score"] == "phosphors_site_probability"

    def test_missing_extra_returns_none(self):
        from qpx.converters.mappings import get_extra

        result = get_extra("diann", "phospho_site_columns")
        assert result is None
