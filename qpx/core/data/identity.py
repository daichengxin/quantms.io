"""Opaque identity hashing for QPX view primary keys.

The ``feature``, ``psm`` and ``pg`` views each carry a single mandatory
identity column (``feature_id`` / ``psm_id`` / ``pg_id``) that is the primary
key of the view. The id is an **opaque** signed 64-bit integer derived
deterministically by the *writer* from a footer-declared composite of existing
columns (the view's ``identity_composite``). It is not meant to be parsed or
reversed — it is only an identity/equality token.

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
    """Recursively normalize a composite value into a JSON-serializable, order-
    insensitive form. List/tuple values (e.g. ``grouped_runs``, ``scan``) are
    **sorted** — by each element's canonical JSON — so ``["r1", "r2"]`` and
    ``["r2", "r1"]`` normalize identically, preserving the order-insensitive
    identity the composite primary key had before it became a single hashed id.
    """
    if isinstance(value, (list, tuple)):
        return sorted((_normalize(v) for v in value), key=lambda x: json.dumps(x, sort_keys=True))
    return value


def canonical(values: list) -> bytes:
    """Encode composite *values* into a single, injective byte string.

    Uses canonical JSON (fixed separators, sorted keys). JSON quotes and escapes
    strings, so — unlike a plain delimiter join — distinct composites can never
    alias to the same bytes even when list elements or fields contain the
    delimiters themselves (e.g. run names with commas/brackets): ``["a,b", "c"]``
    and ``["a", "b,c"]`` encode differently. ``None`` becomes JSON ``null``.
    """
    return json.dumps([_normalize(v) for v in values], separators=(",", ":"), sort_keys=True, ensure_ascii=True).encode("utf-8")


def derive_id(values: list) -> int:
    """Derive a signed 64-bit opaque identity from composite *values*.

    Uses BLAKE2b truncated to 8 bytes. The identity is opaque; uniqueness is
    guaranteed not by construction but by the primary-key uniqueness validation.
    """
    return int.from_bytes(hashlib.blake2b(canonical(values), digest_size=8).digest(), "big", signed=True)
