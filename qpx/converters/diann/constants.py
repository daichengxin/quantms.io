"""DIA-NN converter constants — PTM parsing."""

from __future__ import annotations

import re
from functools import lru_cache

from qpx.converters.ptm import build_proforma, from_proforma

# ---------------------------------------------------------------------------
# PTM parsing: DIA-NN Modified.Sequence -> ProForma
# ---------------------------------------------------------------------------


@lru_cache(maxsize=8192)
def to_proforma(modified_sequence: str) -> str:
    """Convert a DIA-NN Modified.Sequence to ProForma notation.

    DIA-NN encodes modifications with parentheses::

        PEPTM(UniMod:35)IDE       -> PEPTM[UNIMOD:35]IDE
        (UniMod:1)PEPTIDE         -> [UNIMOD:1]-PEPTIDE

    Args:
        modified_sequence: The ``Modified.Sequence`` value from DIA-NN output.

    Returns:
        ProForma string, e.g. ``PEPTM[UNIMOD:35]IDEK``.
    """
    if not modified_sequence:
        return ""

    mods: list[tuple[int, str]] = []
    plain_chars: list[str] = []
    i = 0
    n = len(modified_sequence)
    seq_pos = 0

    while i < n:
        if modified_sequence[i] == "(":
            end = modified_sequence.index(")", i)
            inner = modified_sequence[i + 1 : end]
            normalised = re.sub(r"(?i)unimod:", "UNIMOD:", inner)
            position = 0 if seq_pos == 0 else seq_pos
            mods.append((position, normalised))
            i = end + 1
        else:
            plain_chars.append(modified_sequence[i])
            seq_pos += 1
            i += 1

    plain_seq = "".join(plain_chars)
    return build_proforma(plain_seq, mods)


def to_modifications(modified_sequence: str, sequence: str) -> list[dict] | None:
    """Parse modifications from a DIA-NN Modified.Sequence string.

    Converts to ProForma first, then delegates to the shared ``from_proforma``
    parser to produce the QPX modification list.

    Args:
        modified_sequence: The ``Modified.Sequence`` value from DIA-NN output.
        sequence: Stripped peptide sequence (no modification annotations).

    Returns:
        List of modification dicts per QPX schema, or ``None`` if unmodified.
    """
    proforma = to_proforma(modified_sequence)
    return from_proforma(proforma, sequence, meta=None)
