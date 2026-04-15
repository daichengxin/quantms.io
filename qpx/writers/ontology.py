"""OntologyWriter — writes ontology.parquet files."""

from qpx.core.data import OntologySchema
from qpx.writers.base import BaseWriter


class OntologyWriter(BaseWriter):
    _schema_class = OntologySchema
