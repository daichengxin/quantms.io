"""Tests for field_ontology_entries() with source provenance."""


def test_field_ontology_entries():
    """field_ontology_entries: with source provenance, missing CV, backward compat."""
    from qpx.core.scores import field_ontology_entries

    # With source provenance
    resolved = {"intensity": "Precursor.Quantity", "rt": "RT"}
    entries = field_ontology_entries(
        view="feature",
        resolved_mappings=resolved,
        tool_name="DIA-NN",
    )
    rt_entries = [e for e in entries if e["field_name"] == "rt"]
    assert len(rt_entries) == 1
    assert rt_entries[0]["source_column_name"] == "RT"
    assert rt_entries[0]["source_tool"] == "DIA-NN"
    assert rt_entries[0]["ontology_accession"] == "MS:1000016"
    assert rt_entries[0]["view"] == "feature"

    # Missing CV still written
    resolved = {"lfq": "Precursor.Normalised"}
    entries = field_ontology_entries(
        view="feature",
        resolved_mappings=resolved,
        tool_name="DIA-NN",
    )
    lfq_entries = [e for e in entries if e["field_name"] == "lfq"]
    assert len(lfq_entries) == 1
    assert lfq_entries[0]["source_column_name"] == "Precursor.Normalised"
    assert lfq_entries[0]["ontology_accession"] is None

    # Backward compat (no resolved_mappings)
    entries = field_ontology_entries(view="psm")
    field_names = {e["field_name"] for e in entries}
    assert "posterior_error_probability" in field_names
    assert "rt" in field_names
    for e in entries:
        assert "source_column_name" in e
        assert "source_tool" in e
