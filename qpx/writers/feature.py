"""FeatureWriter — writes feature.parquet files."""

import logging
from collections import Counter
from collections.abc import Sequence

from qpx.core.data import FeatureSchema
from qpx.core.data.identity import derive_id
from qpx.writers.base import BaseWriter

_log = logging.getLogger(__name__)


class FeatureWriter(BaseWriter):
    _schema_class = FeatureSchema

    def __init__(self, *args, override_provided_ids: bool = True, **kwargs):
        """Derive identified Feature ids by default from the declared composite."""
        super().__init__(*args, override_provided_ids=override_provided_ids, **kwargs)
        self._unidentified_fallback_count = 0
        # ids derived for unidentified rows via the composite fallback (no producer
        # id), plus a count of every derived id — used to attribute a duplicate-PK
        # error to the fallback path only when a fallback id is the duplicate one,
        # not when the collision is among identified Features.
        self._fallback_ids: set[int] = set()
        self._all_id_counts: Counter = Counter()

    def _should_override_provided_id(self, index: int, composite_values: dict[str, list]) -> bool:
        """Use the natural composite only for identified Features."""
        return self._override_provided_ids and bool(composite_values["peptidoform"][index])

    def _fallback_collision_error(self, exc: ValueError) -> ValueError:
        """Add remediation only when the duplicate primary key is a fallback-derived id.

        A duplicate among *identified* Features must keep the original error — the
        fallback remediation would wrongly imply the composite fallback is at fault.
        """
        message = str(exc)
        if "primary key" in message.lower() and "duplicate" in message.lower() and self._fallback_derived_duplicate():
            return ValueError(
                f"{message}. Unidentified Feature composite fallback was used for "
                f"{self._unidentified_fallback_count} row(s); provide producer-supplied feature_id values "
                "or use a converter with a reliable identity source"
            )
        return exc

    def _fallback_derived_duplicate(self) -> bool:
        """Whether any fallback-derived id was itself derived for more than one row."""
        return any(self._all_id_counts[fallback_id] > 1 for fallback_id in self._fallback_ids)

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
        self._unidentified_fallback_count += sum(
            not peptidoform and existing[index] is None for index, peptidoform in enumerate(peptidoforms)
        )
        ids = super()._identity_values(existing, composite_values, composite, cv_lists, id_field)
        self._all_id_counts.update(i for i in ids if i is not None)
        for index, peptidoform in enumerate(peptidoforms):
            if not peptidoform and existing[index] is None:
                self._fallback_ids.add(ids[index])
        if not self._override_provided_ids:
            return ids

        for index, peptidoform in enumerate(peptidoforms):
            if peptidoform:
                continue
            provided = existing[index]
            if provided is None:
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

    def write_table(self, table) -> None:
        """Write a table and explain any fallback-derived identity collision."""
        try:
            super().write_table(table)
        except ValueError as exc:
            enriched = self._fallback_collision_error(exc)
            if enriched is exc:
                raise
            raise enriched from exc

    def _validate_identity_uniqueness(self) -> None:
        """Accept unidentified composite fallback only after full-file validation."""
        try:
            super()._validate_identity_uniqueness()
        except ValueError as exc:
            enriched = self._fallback_collision_error(exc)
            if enriched is exc:
                raise
            raise enriched from exc

        if self._unidentified_fallback_count:
            composite = ", ".join(self._identity_composite or ())
            _log.warning(
                "%d unidentified Feature row(s) without a producer feature_id used identity_composite [%s]; "
                "full-file primary-key uniqueness was verified",
                self._unidentified_fallback_count,
                composite,
            )
