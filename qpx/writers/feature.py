"""FeatureWriter — writes feature.parquet files."""

from qpx.core.data import FeatureSchema
from qpx.writers.base import BaseWriter


class FeatureWriter(BaseWriter):
    _schema_class = FeatureSchema
