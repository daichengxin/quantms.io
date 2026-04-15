"""DatasetWriter — writes dataset.parquet files."""

from qpx.core.data import DatasetSchema
from qpx.writers.base import BaseWriter


class DatasetWriter(BaseWriter):
    _schema_class = DatasetSchema
