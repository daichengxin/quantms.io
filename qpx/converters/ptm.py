"""Shared PTM utilities: UNIMOD tables, ProForma builder/parser, precursor m/z.

This flat module consolidates shared logic previously spread across
``qpx.converters.ptm.unimod``, ``qpx.converters.ptm.proforma``,
``qpx.converters.ptm.modifications``, and ``qpx.converters.ptm.__init__``.

Tool-specific constants (MaxQuant _MQ_MOD_TO_UNIMOD, _MQ_SHORT_TO_UNIMOD,
resolve_to_unimod) are **not** included here -- they belong in the respective
converter's ``constants.py``.
"""

from __future__ import annotations

import re
from functools import lru_cache
from typing import Optional

# ---------------------------------------------------------------------------
# UNIMOD mass registry
# ---------------------------------------------------------------------------

UNIMOD_MASS: dict[int, float] = {
    1: 42.010565,
    4: 57.021464,
    5: 43.005814,
    7: 0.984016,
    21: 79.966331,
    28: -17.026549,
    34: 14.015650,
    35: 15.994915,
    36: 28.031300,
    37: 42.046950,
    121: 114.042927,
    122: 27.994915,
    214: -18.010565,
    312: 229.162932,
    385: 304.207146,
    730: 304.207146,
    737: 229.162932,
}

_MASS_TOL = 0.01

_MASS_TO_UNIMOD: list[tuple[float, int]] = sorted((mass, uid) for uid, mass in UNIMOD_MASS.items())


def mass_to_unimod(mass: float) -> Optional[int]:
    """Find the UNIMOD accession whose delta is closest to *mass*.

    Returns ``None`` if no match within ``_MASS_TOL``.
    """
    best_uid = None
    best_diff = _MASS_TOL + 1
    for ref_mass, uid in _MASS_TO_UNIMOD:
        diff = abs(ref_mass - mass)
        if diff < best_diff:
            best_diff = diff
            best_uid = uid
    return best_uid if best_diff <= _MASS_TOL else None


# Monoisotopic amino acid residue masses (Da)
AA_MASS: dict[str, float] = {
    "A": 71.03711,
    "R": 156.10111,
    "N": 114.04293,
    "D": 115.02694,
    "C": 103.00919,
    "E": 129.04259,
    "Q": 128.05858,
    "G": 57.02146,
    "H": 137.05891,
    "I": 113.08406,
    "L": 113.08406,
    "K": 128.09496,
    "M": 131.04049,
    "F": 147.06841,
    "P": 97.05276,
    "S": 87.03203,
    "T": 101.04768,
    "W": 186.07931,
    "Y": 163.06333,
    "V": 99.06841,
}

WATER_MASS = 18.010565
PROTON_MASS = 1.007276


def parse_unimod_delta(mod_str: str) -> Optional[float]:
    """Extract the mass delta from a modification token.

    Accepts:
        - ``UniMod:35`` / ``UNIMOD:35`` -> lookup in ``UNIMOD_MASS``
        - ``+15.9949`` / ``-17.0265`` -> parse directly
    """
    m = re.match(r"(?:UniMod|UNIMOD):(\d+)", mod_str, re.IGNORECASE)
    if m:
        uid = int(m.group(1))
        return UNIMOD_MASS.get(uid)
    try:
        return float(mod_str)
    except ValueError:
        return None


# ---------------------------------------------------------------------------
# ProForma builder
# ---------------------------------------------------------------------------


def build_proforma(sequence: str, mods: list[tuple[int, str]]) -> str:
    """Build a ProForma string from sequence and modification list.

    Args:
        sequence: Stripped peptide sequence (no modification annotations).
        mods: List of (position, tag) where position is 0 for N-term,
            1..N for residue positions (1-indexed). Tag is e.g. "UNIMOD:35"
            or "CHEMMOD:+15.99".

    Returns:
        ProForma string, e.g. ``PEPTM[UNIMOD:35]IDEK`` or ``[UNIMOD:1]-PEPTIDEK``.
    """
    if not sequence:
        return ""

    pos_to_tag: dict[int, str] = dict(mods)
    parts: list[str] = []

    if 0 in pos_to_tag:
        parts.append(f"[{pos_to_tag[0]}]-")

    for i, aa in enumerate(sequence, start=1):
        parts.append(aa)
        if i in pos_to_tag:
            parts.append(f"[{pos_to_tag[i]}]")

    return "".join(parts)


# ---------------------------------------------------------------------------
# ProForma parser -> list<modification>
# ---------------------------------------------------------------------------


def _normalize_peptidoform(peptidoform: str) -> str:
    """Normalise mzTab-style ``(ModName)`` to ProForma ``[ModName]``.

    mzTab encodes modifications with parentheses (e.g.
    ``C(Carbamidomethyl)R``), whereas ProForma uses square brackets.
    This function converts parenthetical mod tags that follow an amino-acid
    residue (or appear at position 0 for N-term mods) to bracket notation
    so that :func:`_from_proforma_impl` can parse them.
    """
    if "(" not in peptidoform:
        return peptidoform
    out: list[str] = []
    n = len(peptidoform)
    i = 0
    while i < n:
        if peptidoform[i] == "(":
            # Find matching closing parenthesis respecting nesting depth
            depth = 1
            j = i + 1
            while j < n and depth > 0:
                if peptidoform[j] == "(":
                    depth += 1
                elif peptidoform[j] == ")":
                    depth -= 1
                j += 1
            if depth != 0:
                out.append(peptidoform[i:])
                break
            # Convert outer delimiters; keep inner content as-is
            out.append("[")
            out.append(peptidoform[i + 1 : j - 1])
            out.append("]")
            i = j
        else:
            out.append(peptidoform[i])
            i += 1
    return "".join(out)


def _from_proforma_impl(
    peptidoform: str,
    sequence: str,
    meta: Optional[dict] = None,
    site_scores: Optional[dict[int, list[dict]]] = None,
) -> Optional[list[dict]]:
    """Core implementation of ProForma modification parsing.

    See :func:`from_proforma` for full documentation.
    """
    # Normalise mzTab parenthetical notation to ProForma brackets
    peptidoform = _normalize_peptidoform(peptidoform)
    if not peptidoform or peptidoform == sequence:
        return None

    mods: dict[str, dict] = {}
    seq_pos = 0
    n = len(peptidoform)

    i = 0
    while i < n:
        if peptidoform[i] == "[":
            try:
                end = peptidoform.index("]", i)
            except ValueError:
                return None  # Malformed ProForma
            mod_str = peptidoform[i + 1 : end]

            position = seq_pos
            aa = sequence[seq_pos - 1] if 0 < seq_pos <= len(sequence) else None

            name = mod_str
            accession = None

            if mod_str.startswith("UNIMOD:"):
                accession = mod_str
                if meta:
                    entry = meta.get(accession)
                    if entry:
                        name = entry[0]
                else:
                    name = accession
            elif mod_str.startswith("+") or mod_str.startswith("-"):
                if meta:
                    for acc, (mod_name, sites, _) in meta.items():
                        if aa in sites or "X" in sites:
                            name = mod_name
                            accession = acc
                            break
                    else:
                        name = f"CHEMMOD:{mod_str}"
                else:
                    name = f"CHEMMOD:{mod_str}"

            key = accession or name
            if key not in mods:
                mods[key] = {"name": name, "accession": accession, "positions": []}
            pos_scores = site_scores.get(position) if site_scores else None
            mods[key]["positions"].append({"position": position, "amino_acid": aa, "scores": pos_scores})

            i = end + 1
        elif peptidoform[i] == "-":
            i += 1
        else:
            seq_pos += 1
            i += 1

    return list(mods.values()) if mods else None


@lru_cache(maxsize=8192)
def _from_proforma_cached(peptidoform: str, sequence: str) -> Optional[list[dict]]:
    """Cached fast path for from_proforma when no meta or site_scores."""
    return _from_proforma_impl(peptidoform, sequence)


def from_proforma(
    peptidoform: str,
    sequence: str,
    meta: Optional[dict] = None,
    site_scores: Optional[dict[int, list[dict]]] = None,
) -> Optional[list[dict]]:
    """Parse modifications from a ProForma-style peptidoform string.

    Handles: ``M[UNIMOD:35]PEPTIDEK``, ``M[+15.9949]PEPTIDEK``,
    ``[UNIMOD:1]-PEPTIDEK``

    Uses an LRU cache for the common case (no *meta*, no *site_scores*).
    Callers must **not** mutate the returned list or its nested dicts, as
    the cached path returns the same object for repeated inputs.

    Args:
        peptidoform: ProForma string with [tag] annotations.
        sequence: Stripped (unmodified) sequence.
        meta: Optional mzTab modifications metadata:
            ``{accession: (name, sites, positions)}``.
            When provided, used to resolve UNIMOD accessions to human-readable
            names and to match mass shifts by site. When None, uses the tag as
            name.
        site_scores: Optional per-position site localization scores.
            Dict mapping position (1-indexed, 0 for N-term) to a list of
            score dicts ``[{score_name, score_value, higher_better}]``.
            Used for phospho site localization probabilities.

    Returns:
        List of modification dicts (``{name, accession, positions}``) per QPX
        schema, or ``None`` if no modifications.
    """
    if meta is None and site_scores is None:
        return _from_proforma_cached(peptidoform, sequence)
    return _from_proforma_impl(peptidoform, sequence, meta, site_scores)


# ---------------------------------------------------------------------------
# Precursor m/z computation
# ---------------------------------------------------------------------------


def compute_precursor_mz(modified_seq: str, charge: int) -> Optional[float]:
    """Compute theoretical precursor m/z for a modified peptide.

    Uses monoisotopic amino-acid masses plus modification deltas encoded as
    ``(UniMod:N)`` tokens.

    Args:
        modified_seq: Modified sequence, e.g. ``PEPTM(UniMod:35)IDE``.
        charge: Precursor charge state.

    Returns:
        Theoretical m/z or ``None`` when unknown modification or charge <= 0.
    """
    if not modified_seq or charge <= 0:
        return None

    mass = 0.0
    i = 0
    n = len(modified_seq)

    while i < n:
        if modified_seq[i] == "(":
            end = modified_seq.index(")", i)
            mod_token = modified_seq[i + 1 : end]
            delta = parse_unimod_delta(mod_token)
            if delta is None:
                return None
            mass += delta
            i = end + 1
        else:
            aa = modified_seq[i]
            aa_mass = AA_MASS.get(aa.upper())
            if aa_mass is None:
                return None
            mass += aa_mass
            i += 1

    neutral_mass = mass + WATER_MASS
    return (neutral_mass + charge * PROTON_MASS) / charge


__all__ = [
    "AA_MASS",
    "PROTON_MASS",
    "UNIMOD_MASS",
    "WATER_MASS",
    "build_proforma",
    "compute_precursor_mz",
    "from_proforma",
    "mass_to_unimod",
    "parse_unimod_delta",
]
