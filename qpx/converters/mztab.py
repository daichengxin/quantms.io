"""Shared mzTab loading logic for DuckDB-based converters.

This module provides helpers that parse an mzTab file (metadata, PSM,
protein, and MSstats sections) and load them into DuckDB tables.  The
quantms adapters compose these helpers rather than inheriting from the
legacy ``MzTabIndexer``.
"""

from __future__ import annotations

import gzip
import logging
import os
import re
import tempfile
from pathlib import Path

import duckdb
import pandas as pd

from qpx.core.sql import escape_path, sql_build, validate_identifier, validate_table

logger = logging.getLogger(__name__)

# Prefixes that identify mzTab line types
_METADATA_PREFIX = "MTD"
_PROTEIN_HEADER_PREFIX = "PRH"
_PROTEIN_LINE_PREFIX = "PRT"
_PSM_HEADER_PREFIX = "PSH"
_PSM_LINE_PREFIX = "PSM"

# Files larger than this threshold use the fast DuckDB-native path
_FAST_LOAD_THRESHOLD_BYTES = 500 * 1024 * 1024  # 500 MB


def _open_mztab(path: str):
    """Open an mzTab file (plain or gzipped)."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


# ------------------------------------------------------------------
# Public helpers
# ------------------------------------------------------------------


def load_mztab_sections(
    conn: duckdb.DuckDBPyConnection,
    mztab_path: str,
) -> None:
    """Parse an mzTab file and load metadata, proteins, and PSMs into DuckDB.

    For large files (>500 MB), uses a fast path that splits the mzTab into
    temporary section files and loads them via DuckDB's native ``read_csv``,
    bypassing the Python → pandas → DuckDB bottleneck.

    After calling this function the connection will contain:
        * ``metadata``  -- two-column table (key TEXT, value TEXT)
        * ``proteins``  -- protein section with dynamic columns
        * ``psms``      -- PSM section with dynamic columns

    Args:
        conn: An open DuckDB connection (in-memory or persistent).
        mztab_path: Path to the mzTab file (plain-text or ``.gz``).
    """
    file_size = os.path.getsize(mztab_path)
    is_gz = str(mztab_path).endswith(".gz")

    if not is_gz and file_size >= _FAST_LOAD_THRESHOLD_BYTES:
        logger.info("Using fast DuckDB-native loader for large mzTab (%.1f GB)", file_size / (1024**3))
        _load_mztab_fast(conn, mztab_path)
    else:
        _load_mztab_classic(conn, mztab_path)


def _load_mztab_classic(
    conn: duckdb.DuckDBPyConnection,
    mztab_path: str,
) -> None:
    """Original in-memory mzTab loader (good for files < 500 MB)."""
    metadata_rows: list[tuple[str, str]] = []
    protein_header: list[str] | None = None
    protein_rows: list[list[str]] = []
    psm_header: list[str] | None = None
    psm_rows: list[list[str]] = []

    with _open_mztab(mztab_path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue

            parts = line.split("\t")
            prefix = parts[0] if parts else ""

            if prefix == _METADATA_PREFIX and len(parts) >= 3:
                metadata_rows.append((parts[1], parts[2]))

            elif prefix == _PROTEIN_HEADER_PREFIX:
                protein_header = [_clean_col(c) for c in parts[1:]]

            elif prefix == _PROTEIN_LINE_PREFIX and protein_header is not None:
                protein_rows.append(parts[1:])

            elif prefix == _PSM_HEADER_PREFIX:
                psm_header = [_clean_col(c) for c in parts[1:]]

            elif prefix == _PSM_LINE_PREFIX and psm_header is not None:
                psm_rows.append(parts[1:])

    _register_metadata(conn, metadata_rows)
    _register_section_df(conn, "proteins", protein_header, protein_rows, "accession")
    _register_section_df(conn, "psms", psm_header, psm_rows, "sequence")

    logger.info(
        "mzTab loaded: %d metadata rows, %d proteins, %d PSMs",
        len(metadata_rows),
        len(protein_rows),
        len(psm_rows),
    )


def _load_mztab_fast(
    conn: duckdb.DuckDBPyConnection,
    mztab_path: str,
) -> None:
    """Fast mzTab loader that splits sections to temp files and uses DuckDB read_csv.

    Single-pass over the mzTab:
      - Metadata rows are small, accumulated in memory.
      - PSM and protein rows are streamed to temporary TSV files.
      - DuckDB reads the temp files natively (no Python/pandas overhead).
    """
    metadata_rows: list[tuple[str, str]] = []
    psm_count = 0
    protein_count = 0

    tmpdir = tempfile.mkdtemp(prefix="mztab_split_")
    psm_tmp = os.path.join(tmpdir, "psms.tsv")
    prot_tmp = os.path.join(tmpdir, "proteins.tsv")

    psm_header_line: str | None = None
    prot_header_line: str | None = None

    try:
        with (
            _open_mztab(mztab_path) as fh,
            open(psm_tmp, "w", encoding="utf-8") as psm_out,
            open(prot_tmp, "w", encoding="utf-8") as prot_out,
        ):
            for line in fh:
                line = line.rstrip("\n\r")
                if not line:
                    continue

                # Check only the first 3 chars for speed
                prefix = line[:3]

                if prefix == _METADATA_PREFIX:
                    parts = line.split("\t", 3)
                    if len(parts) >= 3:
                        metadata_rows.append((parts[1], parts[2]))

                elif prefix == _PSM_HEADER_PREFIX:
                    # Write cleaned header to PSM temp file
                    parts = line.split("\t")
                    cleaned = [_clean_col(c) for c in parts[1:]]
                    psm_header_line = "\t".join(cleaned) + "\n"
                    psm_out.write(psm_header_line)

                elif prefix == _PSM_LINE_PREFIX and psm_header_line is not None:
                    # Strip the "PSM\t" prefix and write the rest directly
                    psm_out.write(line[4:] + "\n")
                    psm_count += 1

                elif prefix == _PROTEIN_HEADER_PREFIX:
                    parts = line.split("\t")
                    cleaned = [_clean_col(c) for c in parts[1:]]
                    prot_header_line = "\t".join(cleaned) + "\n"
                    prot_out.write(prot_header_line)

                elif prefix == _PROTEIN_LINE_PREFIX and prot_header_line is not None:
                    prot_out.write(line[4:] + "\n")
                    protein_count += 1

        _register_metadata(conn, metadata_rows)
        _register_section_csv(conn, "proteins", prot_tmp, prot_header_line is not None, protein_count, "accession")
        _register_section_csv(conn, "psms", psm_tmp, psm_header_line is not None, psm_count, "sequence")

        logger.info(
            "mzTab loaded (fast): %d metadata rows, %d proteins, %d PSMs",
            len(metadata_rows),
            protein_count,
            psm_count,
        )
    finally:
        # Clean up temp files
        for f in (psm_tmp, prot_tmp):
            try:
                os.unlink(f)
            except OSError:
                pass
        try:
            os.rmdir(tmpdir)
        except OSError:
            pass


def load_msstats(
    conn: duckdb.DuckDBPyConnection,
    msstats_path: str,
) -> None:
    """Load an MSstats input CSV/TSV into DuckDB as the ``msstats`` table.

    Args:
        conn: An open DuckDB connection that already has ``metadata``.
        msstats_path: Path to the MSstats input file.
    """
    conn.execute(
        sql_build(
            """CREATE OR REPLACE TABLE msstats AS
        SELECT * FROM read_csv_auto('$path',
            header=true, auto_detect=true,
            quote='"', delim=',', null_padding=true,
            max_line_size=10000000)""",
            path=escape_path(msstats_path),
        )
    )
    count = conn.execute("SELECT COUNT(*) FROM msstats").fetchone()[0]
    logger.info(f"Loaded {count:,} MSstats rows from {msstats_path}")


def extract_ms_runs(conn: duckdb.DuckDBPyConnection) -> dict[int, str]:
    """Return a mapping of MS run index to file stem from metadata.

    Example: ``{1: 'sample_01', 2: 'sample_02'}``
    """
    rows = conn.execute("""
        SELECT key, value FROM metadata
        WHERE key LIKE 'ms_run[%]-location'
        """).fetchall()

    ms_runs: dict[int, str] = {}
    for key, value in rows:
        m = re.search(r"\[(\d+)\]", key)
        if m:
            idx = int(m.group(1))
            ms_runs[idx] = Path(value).stem
    return ms_runs


def extract_modifications(conn: duckdb.DuckDBPyConnection) -> dict:
    """Parse modification metadata from the mzTab metadata table.

    Returns a dict keyed by modification accession, with values
    ``(name, list[site], list[position])``.
    """
    rows = conn.execute("""
        SELECT key, value FROM metadata
        WHERE key LIKE '%_mod[%'
           OR key LIKE '%variable_mod[%'
        ORDER BY key
        """).fetchall()

    temp: dict[str, dict] = {}
    for key, value in rows:
        base_key = key.split("-")[0]
        if base_key not in temp:
            temp[base_key] = {
                "name": None,
                "accession": None,
                "site": None,
                "position": None,
            }

        if "-site" in key:
            temp[base_key]["site"] = value.strip()
        elif "-position" in key:
            temp[base_key]["position"] = value.strip()
        else:
            # Main modification line -- parse [CV, accession, name, ...]
            vals = value.replace("[", "").replace("]", "").split(",")
            if len(vals) >= 3:
                temp[base_key]["accession"] = vals[1].strip()
                temp[base_key]["name"] = vals[2].strip()

    # Group by accession
    modifications: dict = {}
    for entry in temp.values():
        acc = entry["accession"]
        if acc and acc not in modifications:
            modifications[acc] = (
                entry["name"],
                [entry["site"]],
                [entry["position"]],
            )
        elif acc:
            modifications[acc][1].append(entry["site"])
            modifications[acc][2].append(entry["position"])

    return modifications


def extract_score_names(
    conn: duckdb.DuckDBPyConnection,
) -> dict[str, dict[int, str]]:
    """Parse search engine score metadata.

    Returns::

        {"psms": {1: "percolator:score", ...},
         "peptides": {...},
         "proteins": {...}}
    """
    rows = conn.execute("""
        SELECT key, value FROM metadata
        WHERE key LIKE '%search_engine_score%'
        """).fetchall()

    result: dict[str, dict[int, str]] = {
        "psms": {},
        "peptides": {},
        "proteins": {},
    }
    for key, value in rows:
        m = re.search(r"\[(\d+)\]", key)
        if not m:
            continue
        idx = int(m.group(1))
        # Extract the name term (third comma-separated value)
        parts = value.replace("[", "").replace("]", "").split(",")
        name = parts[2].strip() if len(parts) >= 3 else value.strip()

        if key.startswith("psm_"):
            result["psms"][idx] = name
        elif key.startswith("peptide_"):
            result["peptides"][idx] = name
        elif key.startswith("protein_"):
            result["proteins"][idx] = name

    return result


# ------------------------------------------------------------------
# Private helpers
# ------------------------------------------------------------------


def _clean_col(name: str) -> str:
    """Normalise a column name from the mzTab header."""
    return name.strip().lower().replace(" ", "_")


def _register_metadata(
    conn: duckdb.DuckDBPyConnection,
    metadata_rows: list[tuple[str, str]],
) -> None:
    """Load metadata rows into the ``metadata`` DuckDB table."""
    meta_df = pd.DataFrame(metadata_rows, columns=["key", "value"])
    conn.execute("DROP TABLE IF EXISTS metadata")
    conn.from_df(meta_df).create("metadata")


def _register_section_df(
    conn: duckdb.DuckDBPyConnection,
    table_name: str,
    header: list[str] | None,
    rows: list[list[str]],
    fallback_col: str,
) -> None:
    """Load a parsed section into DuckDB, or create an empty fallback table."""
    tbl = validate_table(table_name)
    if header and rows:
        df = pd.DataFrame(rows, columns=header)
        conn.execute(sql_build("DROP TABLE IF EXISTS $t", t=tbl))
        conn.from_df(df).create(table_name)
    else:
        col = validate_identifier(fallback_col)
        conn.execute(sql_build("CREATE TABLE IF NOT EXISTS $t ($col VARCHAR)", t=tbl, col=col))


def _register_section_csv(
    conn: duckdb.DuckDBPyConnection,
    table_name: str,
    tmp_path: str,
    has_header: bool,
    row_count: int,
    fallback_col: str,
) -> None:
    """Load a temp TSV file into DuckDB, or create an empty fallback table."""
    tbl = validate_table(table_name)
    if has_header and row_count > 0:
        conn.execute(
            sql_build(
                """CREATE OR REPLACE TABLE $tbl AS
            SELECT * FROM read_csv('$path',
                header=true, delim='\t', auto_detect=true,
                all_varchar=true, null_padding=true,
                max_line_size=10000000)""",
                tbl=tbl,
                path=escape_path(tmp_path),
            )
        )
    else:
        conn.execute(
            sql_build(
                "CREATE TABLE IF NOT EXISTS $tbl ($col VARCHAR)",
                tbl=tbl,
                col=validate_identifier(fallback_col),
            )
        )
