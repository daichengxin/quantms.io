"""Shared helpers for QPX dataset transforms."""

from pathlib import Path

_QPX_PARQUET_SUFFIXES = (
    ".psm.parquet",
    ".feature.parquet",
    ".pg.parquet",
    ".sample.parquet",
    ".run.parquet",
    ".dataset.parquet",
    ".ontology.parquet",
    ".provenance.parquet",
    ".pepmap.parquet",
    ".mz.parquet",
)


def discover_qpx_file_prefix(dataset_path: Path) -> str:
    """Return the single file prefix shared by QPX Parquet views in a directory."""
    prefixes = {
        path.name[: -len(suffix)]
        for suffix in _QPX_PARQUET_SUFFIXES
        for path in dataset_path.glob(f"*{suffix}")
        if path.is_file()
    }
    if not prefixes:
        raise ValueError(f"No QPX Parquet files found in {dataset_path}")
    if len(prefixes) > 1:
        raise ValueError(f"Multiple QPX file prefixes found in {dataset_path}: {', '.join(sorted(prefixes))}")
    return prefixes.pop()
