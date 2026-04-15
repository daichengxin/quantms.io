"""ProvenanceWriter — writes provenance.parquet files."""

from qpx.core.data import ProvenanceSchema
from qpx.writers.base import BaseWriter


class ProvenanceWriter(BaseWriter):
    _schema_class = ProvenanceSchema
