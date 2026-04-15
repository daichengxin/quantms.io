"""PsmWriter — writes psm.parquet files."""

from qpx.core.data import PsmSchema
from qpx.writers.base import BaseWriter


class PsmWriter(BaseWriter):
    _schema_class = PsmSchema
