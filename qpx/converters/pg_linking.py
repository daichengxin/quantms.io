"""Shared feature->pg cross-reference linker (feature.pg_ids).

Every converter that emits both the ``feature`` and ``pg`` views can stamp each
feature's ``pg_ids`` — the direct references to the pg rows that feature belongs
to (its ``pg_id`` values, ``list<int64>``). A feature maps to the pg row(s) whose
``pg_accessions`` membership equals the feature's group AND whose ``grouped_runs``
quantification unit contains the feature's run — one row per label (so a TMT
feature references one pg_id per channel).

The pg_id is derived from the pg identity composite (``pg.yaml``:
``[pg_accessions, grouped_runs, label]``) with the SAME ``derive_id`` +
order-independent hashing the PgWriter applies in ``qpx/writers/base.py``, so the
ids produced here are byte-identical to the ids stamped on the written pg rows
(bigbio/qpx#266). This module centralizes that rule so every converter
(openms-consensus, DIA-NN, MaxQuant, ...) stays byte-identical to the writer and
to each other, and never re-implements the composite/unordered-index logic.
"""

from __future__ import annotations

from qpx.core.data import PgSchema
from qpx.core.data.identity import derive_id

# The pg view's identity composite (pg.yaml: [pg_accessions, grouped_runs, label])
# and the subset of it hashed order-independently. Mirrors the rule the PgWriter
# applies in qpx/writers/base.py (grouped_runs and pg_accessions are membership
# sets) so a pg_id derived here is byte-identical to the id the writer stamps.
_PG_IDENTITY_COMPOSITE = tuple(PgSchema.identity_composite)
_PG_UNORDERED_INDICES = tuple(
    index for index, field in enumerate(_PG_IDENTITY_COMPOSITE) if field in ("grouped_runs", "pg_accessions")
)


def _accession_of(item) -> str | None:
    """Pull the accession string from a pg_accessions element.

    Handles both the pg view's list-of-strings form and the feature view's
    list-of-structs form (``{"accession": ..., "start": ..., ...}`` dicts, or
    objects exposing an ``accession`` attribute).
    """
    if item is None:
        return None
    if isinstance(item, str):
        return item
    if isinstance(item, dict):
        acc = item.get("accession")
    else:
        acc = getattr(item, "accession", None)
    return None if acc is None else str(acc)


def canonical_membership_key(pg_accessions) -> tuple[str, ...]:
    """Order-independent membership key (sorted accession strings) for a group.

    Accepts either a pg record's ``pg_accessions`` (list of strings) or a feature
    record's ``pg_accessions`` (list of ``{accession, ...}`` structs), so the same
    group keyed from either side collides on the same key. Returns ``()`` for a
    null/empty membership.
    """
    if not pg_accessions:
        return ()
    accessions = {acc for acc in (_accession_of(item) for item in pg_accessions) if acc is not None}
    return tuple(sorted(accessions))


def derive_pg_id(pg_accessions, grouped_runs, label) -> int:
    """Derive the pg_id for one pg row's identity composite.

    ``pg_accessions`` must be the pg view's list-of-strings membership and
    ``grouped_runs`` the row's quantification unit — the SAME values the PgWriter
    persists — so the id is byte-identical to the written ``pg_id``.
    """
    composite = {
        "pg_accessions": list(pg_accessions or []),
        "grouped_runs": list(grouped_runs or []),
        "label": label,
    }
    values = [composite[field] for field in _PG_IDENTITY_COMPOSITE]
    return derive_id(values, unordered_list_indices=_PG_UNORDERED_INDICES)


def build_pg_id_lookup(pg_records) -> dict[tuple[tuple[str, ...], str], list[int]]:
    """Build ``(canonical membership key, run) -> [pg_id, ...]`` from pg records.

    The records are first flattened with the SAME per-label explosion the
    ``PgWriter`` applies on write (``intensities`` list -> one row per label; a
    flat scalar-``label`` record passes through), so each derived ``pg_id`` matches
    a written pg row byte-for-byte. Each flat row's ``pg_id`` is derived from its
    identity composite and its ``grouped_runs`` exploded, so a feature in run ``R``
    whose group membership is ``M`` looks up every pg row (one per label) it
    belongs to via ``(M, R)``. Ids are collected in row order and de-duplicated.
    """
    # Imported lazily to avoid importing the writer stack at module load.
    from qpx.writers.pg import _explode_pg_records

    lookup: dict[tuple[tuple[str, ...], str], list[int]] = {}
    for record in _explode_pg_records(list(pg_records)):
        pg_accessions = record.get("pg_accessions")
        grouped_runs = record.get("grouped_runs") or []
        pg_id = derive_pg_id(pg_accessions, grouped_runs, record.get("label"))
        key = canonical_membership_key(pg_accessions)
        for run in grouped_runs:
            ids = lookup.setdefault((key, run), [])
            if pg_id not in ids:
                ids.append(pg_id)
    return lookup


def update_pg_id_lookup(lookup, pg_records) -> dict[tuple[tuple[str, ...], str], list[int]]:
    """Merge a batch of pg records into an existing lookup (streaming converters).

    Lets a batch-oriented adapter accumulate the full feature->pg lookup one batch
    at a time without holding every pg record in memory. De-duplicates ids per key.
    """
    for key, ids in build_pg_id_lookup(pg_records).items():
        existing = lookup.setdefault(key, [])
        for pg_id in ids:
            if pg_id not in existing:
                existing.append(pg_id)
    return lookup


def pg_ids_for_feature(lookup, feature_pg_accessions, run) -> list[int] | None:
    """Look up a feature's ``pg_ids`` from the membership+run lookup.

    Returns ``None`` when the feature has no group membership or its group has no
    pg row for that run — those features get null ``pg_ids`` (ids are never
    fabricated).
    """
    key = canonical_membership_key(feature_pg_accessions)
    if not key:
        return None
    ids = lookup.get((key, run))
    return list(ids) if ids else None
