"""MzWriter — writes mz.parquet files."""

from qpx.core.data import MzSchema
from qpx.writers.base import BaseWriter


class MzWriter(BaseWriter):
    _schema_class = MzSchema
