"""Ontology data structure — field-to-ontology mappings."""

from qpx.core.data.base import BaseStructure
from qpx.core.data.loader import load_schema
from qpx.core.query import _escape_sql_string

OntologySchema = load_schema("ontology")


class Ontology(BaseStructure):
    """Field-to-ontology mappings for self-describing datasets."""

    _schema_class = OntologySchema

    def by_view(self, view_name: str) -> "Ontology":
        """Filter ontology mappings by QPX view name."""
        return self.filter(f"view = '{_escape_sql_string(view_name)}'")
