"""PepMapWriter — writes pepmap.parquet files."""

from qpx.core.data import PepMapSchema
from qpx.writers.base import BaseWriter


class PepMapWriter(BaseWriter):
    _schema_class = PepMapSchema
