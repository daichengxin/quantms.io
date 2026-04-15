"""DatasetMeta data structure — project-level metadata."""

from qpx.core.data.base import BaseStructure
from qpx.core.data.loader import load_schema

DatasetSchema = load_schema("dataset")


class DatasetMeta(BaseStructure):
    """Project-level metadata (single-row structure)."""

    _schema_class = DatasetSchema
