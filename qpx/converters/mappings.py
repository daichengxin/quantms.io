"""Central column-mapping registry, loaded from column_mappings.yaml."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

_YAML_PATH = Path(__file__).parent / "column_mappings.yaml"
_REGISTRY: dict | None = None


def _load() -> dict:
    global _REGISTRY
    if _REGISTRY is None:
        with open(_YAML_PATH, encoding="utf-8") as f:
            _REGISTRY = yaml.safe_load(f)
    return _REGISTRY


def get_field_mappings(tool: str, view: str) -> dict[str, list[str]]:
    """Return column mappings for a tool + view (feature/psm/pg)."""
    registry = _load()
    if tool not in registry:
        raise KeyError(f"Unknown tool '{tool}'. Available: {sorted(registry)}")
    views = registry[tool].get("column_mapping", {})
    if view not in views:
        raise KeyError(f"Tool '{tool}' has no '{view}' view. Available: {sorted(views)}")
    return views[view]


def get_tool_meta(tool: str) -> dict[str, str]:
    """Return tool_name, tool_versions, rt_unit for a tool."""
    entry = _load()[tool]
    return {
        "tool_name": entry["tool_name"],
        "tool_versions": entry["tool_versions"],
        "rt_unit": entry.get("rt_unit", "minute"),
    }


def get_extra(tool: str, key: str) -> Any:
    """Return a tool-specific extra config section, or None if absent."""
    return _load()[tool].get(key)
