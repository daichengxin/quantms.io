"""Provenance data structure — processing pipeline provenance."""

from qpx.core.data.base import BaseStructure
from qpx.core.data.loader import load_schema

ProvenanceSchema = load_schema("provenance")


class Provenance(BaseStructure):
    """Processing provenance — one row per pipeline step."""

    _schema_class = ProvenanceSchema
