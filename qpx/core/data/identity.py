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

# Unit Separator between composite fields, and a NUL sentinel for None, so that
# the canonical byte string is unambiguous (no ordinary field value contains
# these control characters).
FIELD_SEP = "\x1f"
NULL_TOKEN = "\x00"


def canonical(values: list) -> bytes:
    """Encode composite *values* into a single unambiguous byte string."""
    return FIELD_SEP.join(NULL_TOKEN if v is None else str(v) for v in values).encode("utf-8")


def derive_id(values: list) -> int:
    """Derive a signed 64-bit opaque identity from composite *values*.

    Uses BLAKE2b truncated to 8 bytes. The identity is opaque; uniqueness is
    guaranteed not by construction but by the primary-key uniqueness validation.
    """
    return int.from_bytes(hashlib.blake2b(canonical(values), digest_size=8).digest(), "big", signed=True)
