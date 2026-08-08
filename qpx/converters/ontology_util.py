"""Shared helpers for converter ontology writing."""

from __future__ import annotations


def dedupe_ontology_entries(entries: list[dict]) -> list[dict]:
    """Collapse ontology entries to one row per ``(field_name, view)`` primary key.

    Order-preserving, first-wins: later duplicates for the same key are dropped.
    """
    seen: set[tuple] = set()
    deduped: list[dict] = []
    for entry in entries:
        key = (entry.get("field_name"), entry.get("view"))
        if key in seen:
            continue
        seen.add(key)
        deduped.append(entry)
    return deduped
