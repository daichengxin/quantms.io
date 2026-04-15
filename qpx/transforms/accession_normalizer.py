"""Normalize protein accessions between UniProt full format and bare accession.

Forward:  sp|P04114|APOB_HUMAN  →  P04114
Reverse:  P04114  →  sp|P04114|APOB_HUMAN  (requires FASTA database)

Handles prefixes: sp|, tr|, CONTAM_sp|, ENTRAP_sp|
Contaminant/entrapment prefixes (CONTAM_, ENTRAP_) are preserved.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

logger = logging.getLogger(__name__)

# Matches:  sp|ACCESSION|NAME  or  tr|ACCESSION|NAME
# Also:     CONTAM_sp|CONTAM_ACCESSION|...  or  ENTRAP_sp|ENTRAP_ACCESSION|...
_FULL_RE = re.compile(
    r"^(?P<prefix>(?:CONTAM_|ENTRAP_)?(?:sp|tr))\|"
    r"(?P<accession>[^|]+)\|"
    r"(?P<name>.+)$"
)


def to_bare(accession: str) -> str:
    """Convert full UniProt accession to bare form.

    sp|P04114|APOB_HUMAN           → P04114
    CONTAM_sp|CONTAM_P02768|...    → CONTAM_P02768
    tr|A0A123|A0A123_HUMAN         → A0A123
    P04114                         → P04114  (no-op)
    """
    if accession is None:
        return accession
    m = _FULL_RE.match(accession)
    if m:
        return m.group("accession")
    return accession


def to_full(accession: str, lookup: dict[str, str]) -> str:
    """Convert bare accession to full UniProt format using a FASTA-derived lookup.

    P04114 → sp|P04114|APOB_HUMAN  (if in lookup)
    P04114 → P04114                (if not in lookup, returned as-is)
    """
    if accession is None:
        return accession
    # Already in full format?
    if _FULL_RE.match(accession):
        return accession
    return lookup.get(accession, accession)


def parse_fasta_headers(fasta_path: str | Path) -> dict[str, str]:
    """Parse FASTA file and return {bare_accession: full_header} mapping.

    Reads headers like:
        >sp|P04114|APOB_HUMAN Apolipoprotein B-100 ...
        >tr|A0A123|A0A123_HUMAN ...
        >CONTAM_sp|CONTAM_P02768|CONTAM_ALBU_HUMAN ...

    Returns:
        {"P04114": "sp|P04114|APOB_HUMAN",
         "A0A123": "tr|A0A123|A0A123_HUMAN",
         "CONTAM_P02768": "CONTAM_sp|CONTAM_P02768|CONTAM_ALBU_HUMAN"}
    """
    lookup: dict[str, str] = {}
    fasta_path = Path(fasta_path)
    with open(fasta_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            # Extract the identifier (first whitespace-delimited token after >)
            header = line[1:].split()[0]
            m = _FULL_RE.match(header)
            if m:
                lookup[m.group("accession")] = header
    logger.info("Parsed %d protein headers from %s", len(lookup), fasta_path.name)
    return lookup


def normalize_parquet(
    parquet_path: Path,
    output_path: Path,
    direction: str = "forward",
    fasta_lookup: dict[str, str] | None = None,
) -> dict:
    """Normalize accessions in a QPX parquet file.

    Args:
        parquet_path: Input parquet file (feature.parquet or pg.parquet).
        output_path: Where to write the normalized file.
        direction: "forward" (full → bare) or "reverse" (bare → full).
        fasta_lookup: Required for reverse. {bare_accession: full_header}.

    Returns:
        Dict with stats: {"rows", "accessions_changed", "file"}.
    """
    if direction == "reverse" and not fasta_lookup:
        raise ValueError("FASTA lookup required for reverse normalization")

    transform_fn = to_bare if direction == "forward" else lambda a: to_full(a, fasta_lookup)

    table = pq.read_table(str(parquet_path))
    schema = table.schema
    n_changed = 0
    total_rows = len(table)

    # --- Normalize anchor_protein column ---
    if "anchor_protein" in schema.names:
        col = table.column("anchor_protein")
        new_values = []
        for val in col.to_pylist():
            new_val = transform_fn(val) if val else val
            if new_val != val:
                n_changed += 1
            new_values.append(new_val)
        table = table.set_column(
            schema.get_field_index("anchor_protein"),
            "anchor_protein",
            pa.array(new_values, type=pa.string()),
        )

    # --- Normalize pg_accessions column ---
    if "pg_accessions" in schema.names:
        pg_col = table.column("pg_accessions")
        pg_type = pg_col.type

        # Two formats: list<string> (pg.parquet) or list<struct{accession,...}> (feature.parquet)
        if _is_list_of_strings(pg_type):
            new_pg = _normalize_list_of_strings(pg_col, transform_fn)
            table = table.set_column(
                schema.get_field_index("pg_accessions"),
                "pg_accessions",
                new_pg,
            )
        elif _is_list_of_structs_with_accession(pg_type):
            new_pg = _normalize_list_of_structs(pg_col, transform_fn)
            table = table.set_column(
                schema.get_field_index("pg_accessions"),
                "pg_accessions",
                new_pg,
            )

    pq.write_table(table, str(output_path))
    logger.info(
        "Normalized %s → %s (%d rows, %d anchor_protein changes)",
        parquet_path.name,
        output_path.name,
        total_rows,
        n_changed,
    )
    return {"rows": total_rows, "accessions_changed": n_changed, "file": str(output_path)}


def _is_list_of_strings(arrow_type: pa.DataType) -> bool:
    """Check if type is list<string> or large_list<string>."""
    if not (pa.types.is_list(arrow_type) or pa.types.is_large_list(arrow_type)):
        return False
    return pa.types.is_string(arrow_type.value_type) or pa.types.is_large_string(arrow_type.value_type)


def _is_list_of_structs_with_accession(arrow_type: pa.DataType) -> bool:
    """Check if type is list<struct{accession: string, ...}>."""
    if not (pa.types.is_list(arrow_type) or pa.types.is_large_list(arrow_type)):
        return False
    vt = arrow_type.value_type
    if not pa.types.is_struct(vt):
        return False
    return any(f.name == "accession" for f in vt)


def _normalize_list_of_strings(col: pa.ChunkedArray, transform_fn) -> pa.Array:
    """Normalize a list<string> column."""
    new_values = []
    for row in col.to_pylist():
        if row is None:
            new_values.append(None)
        else:
            new_values.append([transform_fn(v) if v else v for v in row])
    return pa.array(new_values, type=col.type)


def _normalize_list_of_structs(col: pa.ChunkedArray, transform_fn) -> pa.Array:
    """Normalize the 'accession' field in a list<struct{accession,...}> column."""
    new_values = []
    for row in col.to_pylist():
        if row is None:
            new_values.append(None)
        else:
            new_row = []
            for item in row:
                if item and "accession" in item:
                    item = dict(item)
                    item["accession"] = transform_fn(item["accession"]) if item["accession"] else item["accession"]
                new_row.append(item)
            new_values.append(new_row)
    return pa.array(new_values, type=col.type)
