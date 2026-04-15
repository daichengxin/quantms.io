"""RunWriter — writes run.parquet files."""

from qpx.core.data import RunSchema
from qpx.writers.base import BaseWriter


class RunWriter(BaseWriter):
    _schema_class = RunSchema
