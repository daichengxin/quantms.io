"""Missed-cleavage computation from peptide sequence + enzyme name.

Provides ``count_missed_cleavages`` which counts internal sites that
match an enzyme's cleavage rule but were not cleaved, and an
``ENZYME_RULES`` lookup table for common proteolytic enzymes.
"""

from __future__ import annotations

import re
from functools import lru_cache

# ---------------------------------------------------------------------------
# Enzyme cleavage rules: compiled regex for internal cleavage sites
# ---------------------------------------------------------------------------

ENZYME_RULES: dict[str, re.Pattern] = {
    # Trypsin: cleaves after K or R, unless followed by P
    "Trypsin": re.compile(r"(?<=[KR])(?!P)"),
    "Trypsin/P": re.compile(r"(?<=[KR])"),
    "Lys-C": re.compile(r"(?<=K)"),
    "Lys-N": re.compile(r"(?=K)"),
    "Arg-C": re.compile(r"(?<=R)"),
    "Asp-N": re.compile(r"(?=[BD])"),
    "Glu-C": re.compile(r"(?<=[DE])"),
    "Chymotrypsin": re.compile(r"(?<=[FYWL])(?!P)"),
}

# Normalised aliases (lower-case) -> canonical name
_ENZYME_ALIASES: dict[str, str] = {}
for _name in ENZYME_RULES:
    _ENZYME_ALIASES[_name.lower()] = _name
# Additional common aliases
_ENZYME_ALIASES["trypsin/p"] = "Trypsin/P"
_ENZYME_ALIASES["lysc"] = "Lys-C"
_ENZYME_ALIASES["lysn"] = "Lys-N"
_ENZYME_ALIASES["argc"] = "Arg-C"
_ENZYME_ALIASES["aspn"] = "Asp-N"
_ENZYME_ALIASES["gluc"] = "Glu-C"


def _resolve_enzyme(enzyme_name: str) -> re.Pattern | None:
    """Resolve an enzyme name (case-insensitive) to its cleavage pattern."""
    canonical = _ENZYME_ALIASES.get(enzyme_name.lower())
    if canonical:
        return ENZYME_RULES[canonical]
    return ENZYME_RULES.get(enzyme_name)


@lru_cache(maxsize=16384)
def count_missed_cleavages(sequence: str, enzyme_name: str) -> int | None:
    """Count missed enzymatic cleavage sites in a peptide sequence.

    A missed cleavage is an internal position in the sequence that matches
    the enzyme's cleavage specificity but was not cleaved (i.e. the peptide
    spans across it).

    Args:
        sequence: Unmodified (stripped) peptide sequence, e.g. ``"PEPTIDEK"``.
        enzyme_name: Enzyme name, e.g. ``"Trypsin"``.  Case-insensitive.

    Returns:
        Number of internal cleavage sites, or ``None`` if the enzyme is
        unknown and cannot be resolved.
    """
    if not sequence:
        return None
    rule = _resolve_enzyme(enzyme_name)
    if rule is None:
        return None
    # Find all positions where the enzyme would cleave
    sites = [m.start() for m in rule.finditer(sequence)]
    # Exclude terminal cleavage (at the very end of the peptide)
    # — a site at position len(sequence) is the C-terminal cut, not missed
    internal = [s for s in sites if 0 < s < len(sequence)]
    return len(internal)
