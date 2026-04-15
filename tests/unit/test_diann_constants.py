"""Tests for DIA-NN column mappings (now loaded from central YAML)."""


def test_diann_mappings():
    """DIA-NN mappings: structure, field mappings are lists, key fields present."""
    from qpx.converters.mappings import get_field_mappings, get_tool_meta

    meta = get_tool_meta("diann")
    assert meta["tool_name"] == "DIA-NN"
    assert isinstance(meta["tool_versions"], str)

    for view in ("feature", "pg"):
        fields = get_field_mappings("diann", view)
        for qpx_field, candidates in fields.items():
            assert isinstance(candidates, list), f"diann.{view}.{qpx_field} must be a list"
            assert len(candidates) > 0

    feature = get_field_mappings("diann", "feature")
    assert "intensity" in feature
    assert "posterior_error_probability" in feature
    assert "rt" in feature
    assert "pg_accessions" in feature
