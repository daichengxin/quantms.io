"""FeatureWriter — writes feature.parquet files."""

from collections.abc import Sequence

from qpx.core.data import FeatureSchema
from qpx.core.data.identity import derive_id
from qpx.writers.base import BaseWriter


class FeatureWriter(BaseWriter):
    _schema_class = FeatureSchema

    def __init__(self, *args, override_provided_ids: bool = True, **kwargs):
        """Derive identified Feature ids by default from the declared composite."""
        super().__init__(*args, override_provided_ids=override_provided_ids, **kwargs)

    def _should_override_provided_id(self, index: int, composite_values: dict[str, list]) -> bool:
        """Use the natural composite only for identified Features."""
        return self._override_provided_ids and bool(composite_values["peptidoform"][index])

    def _identity_values(
        self,
        existing: list,
        composite_values: dict[str, list],
        composite: Sequence[str],
        cv_lists: list | None,
        id_field: str,
    ) -> list[int]:
        """Derive identified ids; for unidentified Features derive from the composite
        when the producer supplies no id, or namespace the producer id when it does."""
        peptidoforms = composite_values["peptidoform"]
        ids = super()._identity_values(existing, composite_values, composite, cv_lists, id_field)
        if not self._override_provided_ids:
            return ids

        for index, peptidoform in enumerate(peptidoforms):
            if peptidoform:
                continue
            provided = existing[index]
            if provided is None:
                # Unidentified Feature with no producer id: the base has already
                # derived an id from the composite (a null peptidoform plus the
                # measured charge/run/rt/…), the best identity available. A residual
                # collision surfaces as a duplicate-primary-key warning at close.
                continue
            # Unidentified Feature WITH a producer id: namespace it by run/unit so
            # tool-local ids cannot collide across runs, and stash the original.
            namespace_field = "run_file_name" if "run_file_name" in composite_values else "quantification_unit_id"
            ids[index] = derive_id([composite_values[namespace_field][index], provided])
            if cv_lists is not None:
                params = list(cv_lists[index] or [])
                provided_param = {"cv_name": f"provided_{id_field}", "cv_value": str(provided)}
                if provided_param not in params:
                    params.append(provided_param)
                cv_lists[index] = params
            self.overridden_id_count += 1
        return ids
