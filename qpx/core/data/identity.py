"""Opaque identity hashing for QPX view primary keys.

The ``feature``, ``psm`` and ``pg`` views each carry a single mandatory
identity column (``feature_id`` / ``psm_id`` / ``pg_id``) that is the primary
key of the view. When no producer ID is supplied, the *writer* derives an
**opaque** signed 64-bit integer deterministically from a footer-declared
composite of existing columns (the view's ``identity_composite``). It is not
meant to be parsed or reversed — it is only an identity/equality token.

Because the id is a fixed-width hash of the composite, distinct composites can
in principle collide onto the same id. That is not silently tolerated: the id
is the primary key, and the primary-key **uniqueness** validation
(``duplicate_pk``) catches any collision, so a real clash surfaces as a
validation issue rather than a silent data loss.
"""

from __future__ import annotations

import hashlib
import json


def _normalize(value):
    """Recursively normalize a composite value into a JSON-serializable form."""
    if isinstance(value, (list, tuple)):
        return [_normalize(v) for v in value]
    return value


def canonical(values: list, *, unordered_list_indices: tuple[int, ...] = ()) -> bytes:
    """Encode composite *values* into a canonical byte string.

    Uses canonical JSON (fixed separators, sorted keys). JSON quotes and escapes
    strings, so — unlike a plain delimiter join — distinct composites can never
    alias to the same bytes even when list elements or fields contain the
    delimiters themselves (e.g. run names with commas/brackets): ``["a,b", "c"]``
    and ``["a", "b,c"]`` encode differently. ``None`` becomes JSON ``null``.
    List order and duplicate elements are preserved unless its top-level
    position is explicitly listed in *unordered_list_indices*. QPX uses that
    exception for set-valued fields such as ``grouped_runs`` and
    ``pg_accessions``; ordered identifiers such as a multi-component ``scan``
    keep their component order.
    """
    normalized = [_normalize(value) for value in values]
    for index in unordered_list_indices:
        value = normalized[index]
        if isinstance(value, list):
            by_encoding = {json.dumps(item, separators=(",", ":"), sort_keys=True, ensure_ascii=True): item for item in value}
            normalized[index] = [by_encoding[key] for key in sorted(by_encoding)]
    return json.dumps(normalized, separators=(",", ":"), sort_keys=True, ensure_ascii=True).encode("utf-8")


def derive_id(values: list, *, unordered_list_indices: tuple[int, ...] = ()) -> int:
    """Derive a signed 64-bit opaque identity from composite *values*.

    Uses BLAKE2b truncated to 8 bytes. The identity is opaque; uniqueness is
    guaranteed not by construction but by the primary-key uniqueness validation.
    """
    encoded = canonical(values, unordered_list_indices=unordered_list_indices)
    return int.from_bytes(hashlib.blake2b(encoded, digest_size=8).digest(), "big", signed=True)
