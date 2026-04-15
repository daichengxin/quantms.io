"""Tests for resolve_columns() utility."""


def test_resolve_columns():
    """resolve_columns: basic, ordered fallback, first match wins, missing skipped, empty, hasattr checks."""
    from qpx.converters.base import BaseConverter, resolve_columns

    # Basic resolution
    mappings = {"intensity": ["Precursor.Quantity"], "rt": ["RT"]}
    available = {"Precursor.Quantity", "RT", "Run", "Sequence"}
    result = resolve_columns(mappings, available)
    assert result == {"intensity": "Precursor.Quantity", "rt": "RT"}

    # Ordered fallback
    result = resolve_columns({"rt": ["RT", "RetentionTime"]}, {"RetentionTime", "Sequence"})
    assert result == {"rt": "RetentionTime"}

    # First match wins
    result = resolve_columns({"rt": ["RT", "RetentionTime"]}, {"RT", "RetentionTime"})
    assert result == {"rt": "RT"}

    # Missing skipped
    result = resolve_columns(
        {"intensity": ["Precursor.Quantity"], "predicted_rt": ["Predicted.RT"]},
        {"Precursor.Quantity"},
    )
    assert result == {"intensity": "Precursor.Quantity"}
    assert "predicted_rt" not in result

    # Empty
    assert resolve_columns({}, {"col1", "col2"}) == {}

    # BaseConverter method existence
    assert hasattr(BaseConverter, "write_ontology")
    assert hasattr(BaseConverter, "write_score_ontology")
