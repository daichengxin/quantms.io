"""File-name helpers shared across QPX inputs."""

from __future__ import annotations

import re
from pathlib import PurePosixPath

_RUN_FILE_SUFFIX = re.compile(r"(?i)\.(?:mzml(?:\.gz)?|mzxml|raw|d|wiff|mgf|dia)$")


def run_file_stem(value: str) -> str:
    """Return a run basename without its acquisition-file suffix."""
    name = PurePosixPath(str(value).strip().replace("\\", "/")).name
    return _RUN_FILE_SUFFIX.sub("", name)
