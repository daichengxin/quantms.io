"""QPX configuration module."""

from pathlib import Path

import yaml

# Default filter configuration path
DEFAULT_FILTERS_PATH = Path(__file__).parent / "default_filters.yaml"


def load_default_filters() -> dict:
    """Load default filter configuration from YAML file."""
    with open(DEFAULT_FILTERS_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


# Lazy loading: only load when accessed
_default_filters = None


def get_default_filters() -> dict:
    """Get default filter configuration (loaded once)."""
    global _default_filters
    if _default_filters is None:
        _default_filters = load_default_filters()
    return _default_filters
