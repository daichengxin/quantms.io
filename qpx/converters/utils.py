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


def cv_column_name(cv_term: str, suffix: str) -> str:
    """Build the mzTab opt_global column name for a CV term.

    Returns the lowercase variant.  Use :func:`cv_column_names` for both
    variants.
    """
    cv_code = cv_term.split(":")[1]
    return f"opt_global_cv_ms:{cv_code}_{suffix}"


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


def clean_peptidoform(peptidoform: str) -> str:
    """Strip leading/trailing underscores from a MaxQuant peptidoform."""
    if not isinstance(peptidoform, str):
        return ""
    return peptidoform.strip("_")
