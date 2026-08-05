"""FeatureWriter — writes feature.parquet files."""

from qpx.core.data import FeatureSchema
from qpx.writers.base import BaseWriter


class FeatureWriter(BaseWriter):
    _schema_class = FeatureSchema

    def __init__(self, *args, override_provided_ids: bool = True, **kwargs):
        """Derive identified Feature ids by default from the declared composite."""
        super().__init__(*args, override_provided_ids=override_provided_ids, **kwargs)

    def _should_override_provided_id(self, index: int, composite_values: dict[str, list]) -> bool:
        """Retain producer ids for unidentified Features."""
        return self._override_provided_ids and bool(composite_values["peptidoform"][index])
