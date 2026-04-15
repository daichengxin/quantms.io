"""PgWriter — writes pg.parquet files."""

from qpx.core.data import PgSchema
from qpx.writers.base import BaseWriter


class PgWriter(BaseWriter):
    _schema_class = PgSchema
