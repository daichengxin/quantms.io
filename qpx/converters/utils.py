from __future__ import annotations

import re
from typing import Optional


def safe_float(val) -> Optional[float]:
    """Convert a value to float, returning None for missing/invalid values.

    Handles None, empty strings, the literal ``"null"``, NaN, and
    anything that cannot be cast to float.
    """
    if val is None or val == "" or val == "null":
        return None
    try:
        f = float(val)
        return f if f == f else None  # NaN check
    except (ValueError, TypeError):
        return None


# ------------------------------------------------------------------
# UniProt accession helpers
# ------------------------------------------------------------------


def parse_uniprot_id(entry: str) -> tuple[str, str]:
    """Split a UniProt-style ``db|ACCESSION|NAME`` id into ``(accession, name)``.

    ``sp|P12345|PROT_HUMAN`` -> ``("P12345", "PROT_HUMAN")``. With only two
    pipe-fields the accession doubles as the name; with none the whole ``entry``
    is both. Shared by the FragPipe/quantms protein-field parsers.
    """
    parts = entry.split("|")
    if len(parts) >= 3:
        return parts[1], parts[2]
    if len(parts) == 2:
        return parts[1], parts[1]
    return entry, entry


def strip_uniprot_prefix(accession: str) -> str:
    """Strip a leading ``sp|``/``tr|`` UniProt db prefix, returning the accession.

    ``sp|P55011|S12A2_HUMAN`` -> ``P55011``; ``tr|A0A..|..._HUMAN`` -> ``A0A..``;
    ids without that prefix (bare accessions, ``CON__``/``REV__`` decoys) are
    returned unchanged. Used by the MaxQuant adapters, whose FASTA prefixing is
    the only case that should be stripped.
    """
    if accession and accession.startswith(("sp|", "tr|")):
        parts = accession.split("|")
        if len(parts) >= 2:
            return parts[1]
    return accession


# ------------------------------------------------------------------
# mzTab spectra_ref helpers
# ------------------------------------------------------------------


def parse_scan_numbers(spectra_ref: str) -> list[int]:
    """Extract scan numbers from an mzTab spectra_ref value.

    Examples:
        ``ms_run[1]:scan=1234``        -> ``[1234]``
        ``ms_run[1]:index=42``         -> ``[42]``
        ``ms_run[1]:spectrum=5``       -> ``[5]``
    """
    if not spectra_ref or spectra_ref == "null":
        return []
    match = re.search(r"(?:scan|index|spectrum)=(\d+)", spectra_ref)
    if match:
        return [int(match.group(1))]
    # Fallback: try extracting any trailing integer after ':'
    parts = spectra_ref.split(":")
    if len(parts) >= 2:
        try:
            return [int(parts[-1])]
        except ValueError:
            pass
    return []


def resolve_run_file(spectra_ref: str, ms_runs: dict[int, str]) -> Optional[str]:
    """Map spectra_ref to its run file stem via ms_runs dict."""
    if not spectra_ref or spectra_ref == "null":
        return None
    m = re.search(r"\[(\d+)\]", spectra_ref)
    if m:
        return ms_runs.get(int(m.group(1)))
    return None


# ------------------------------------------------------------------
# CV term column-name helpers
# ------------------------------------------------------------------


def cv_column_names(cv_term: str, suffix: str) -> tuple[str, str]:
    """Return (lowercase, uppercase) opt_global column names for a CV term."""
    cv_code = cv_term.split(":")[1]
    return (
        f"opt_global_cv_ms:{cv_code}_{suffix}",
        f"opt_global_cv_MS:{cv_code}_{suffix}",
    )


def get_cv_value(row: dict, cv_term: str, suffix: str, default=None):
    """Get a value from a row dict trying both CV column name variants."""
    lo, hi = cv_column_names(cv_term, suffix)
    return row.get(lo, row.get(hi, default))


# ------------------------------------------------------------------
# MaxQuant helpers
# ------------------------------------------------------------------


def mq_flag_to_bool(val) -> bool:
    """Convert MaxQuant '+' flag to boolean."""
    return str(val).strip() == "+"
