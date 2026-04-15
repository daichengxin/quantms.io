"""Sample data structure — biological sample metadata."""

from qpx.core.data.base import BaseStructure
from qpx.core.data.loader import load_schema
from qpx.core.query import _escape_sql_string

SampleSchema = load_schema("sample")


class Sample(BaseStructure):
    """Biological sample metadata."""

    _schema_class = SampleSchema

    def by_organism(self, organism: str) -> "Sample":
        """Filter samples by organism (substring match)."""
        return self.filter(f"organism LIKE '%{_escape_sql_string(organism)}%'")

    def by_disease(self, disease: str) -> "Sample":
        """Filter samples by disease (substring match)."""
        return self.filter(f"disease LIKE '%{_escape_sql_string(disease)}%'")
