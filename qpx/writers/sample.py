"""SampleWriter — writes sample.parquet files."""

from qpx.core.data import SampleSchema
from qpx.writers.base import BaseWriter


class SampleWriter(BaseWriter):
    _schema_class = SampleSchema
