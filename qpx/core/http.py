"""Safe HTTP utilities for URL fetching.

Centralises all ``urlopen`` usage behind scheme-validated helpers so that
Bandit B310 (``audit url open for permitted schemes``) is confined to this
single module.  Callers import :func:`safe_urlopen` instead of using
``urllib.request.urlopen`` directly.
"""

from __future__ import annotations

import urllib.request
from urllib.parse import urlparse

_ALLOWED_SCHEMES = frozenset({"https", "http"})


def safe_urlopen(
    url: str | urllib.request.Request,
    *,
    timeout: int = 30,
    allowed_schemes: frozenset[str] = _ALLOWED_SCHEMES,
):
    """Open *url* after validating its scheme.

    Parameters
    ----------
    url : str or Request
        The URL (or ``Request`` object) to open.
    timeout : int
        Socket timeout in seconds.
    allowed_schemes : frozenset[str]
        Permitted URL schemes (default ``{"https", "http"}``).

    Raises
    ------
    ValueError
        If the URL scheme is not in *allowed_schemes*.
    """
    raw = url.full_url if isinstance(url, urllib.request.Request) else url
    parsed = urlparse(raw)
    if parsed.scheme not in allowed_schemes:
        raise ValueError(f"URL scheme {parsed.scheme!r} not allowed (permitted: {', '.join(sorted(allowed_schemes))}): {raw}")
    return urllib.request.urlopen(url, timeout=timeout)
