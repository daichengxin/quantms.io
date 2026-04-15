"""PRIDE Archive REST API client.

Fetches project metadata from the PRIDE REST API v3
(https://www.ebi.ac.uk/pride/ws/archive/v3/).
"""

from __future__ import annotations

import json
import logging
from urllib.error import HTTPError, URLError
from urllib.request import Request

from qpx.core.http import safe_urlopen

logger = logging.getLogger(__name__)

PRIDE_API_BASE = "https://www.ebi.ac.uk/pride/ws/archive/v3"


def fetch_pride_metadata(accession: str, timeout: int = 30) -> dict:
    """Fetch project metadata from the PRIDE REST API v3.

    Parameters
    ----------
    accession : str
        ProteomeXchange/PRIDE accession (e.g., "PXD054720").
    timeout : int
        HTTP request timeout in seconds.

    Returns
    -------
    dict
        Parsed metadata with keys: project_title, project_description,
        pubmed_id, doi, organisms, instruments, keywords,
        submission_date, publication_date.

    Raises
    ------
    ValueError
        If the accession is not found (HTTP 404).
    ConnectionError
        If the PRIDE API is unreachable.
    """
    url = f"{PRIDE_API_BASE}/projects/{accession}"
    req = Request(url, headers={"Accept": "application/json"})

    try:
        with safe_urlopen(req, timeout=timeout) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except HTTPError as exc:
        if exc.code == 404:
            raise ValueError(f"PRIDE project not found: {accession}") from exc
        raise ConnectionError(f"PRIDE API error (HTTP {exc.code}) for {accession}") from exc
    except URLError as exc:
        raise ConnectionError(f"Cannot reach PRIDE API: {exc.reason}") from exc

    # Extract PubMed ID from references
    pubmed_id = None
    references = data.get("references") or []
    if references:
        pm = references[0].get("pubmedID")
        if pm is not None:
            pubmed_id = str(pm)

    return {
        "project_title": data.get("title"),
        "project_description": data.get("projectDescription"),
        "pubmed_id": pubmed_id,
        "doi": data.get("doi"),
        "organisms": [o.get("name", "") for o in (data.get("organisms") or [])],
        "instruments": [i.get("name", "") for i in (data.get("instruments") or [])],
        "keywords": data.get("keywords") or [],
        "submission_date": data.get("submissionDate"),
        "publication_date": data.get("publicationDate"),
    }
