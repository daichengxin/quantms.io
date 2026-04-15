"""Tests for PRIDE API client."""

import json
from unittest.mock import MagicMock, patch

import pytest

from qpx.core.pride import fetch_pride_metadata

# Sample PRIDE API v3 response
MOCK_PRIDE_RESPONSE = {
    "accession": "PXD014414",
    "title": "Quantitative proteomic landscape of metaplastic breast carcinoma",
    "projectDescription": "This study provides a comprehensive proteomic analysis...",
    "doi": "10.6019/PXD014414",
    "submissionDate": "2019-06-25",
    "publicationDate": "2020-04-08",
    "keywords": ["Proteomics", "Breast cancer"],
    "organisms": [{"@type": "CvParam", "name": "Homo sapiens (human)", "accession": "NEWT:9606"}],
    "instruments": [{"@type": "CvParam", "name": "Q Exactive HF", "accession": "MS:1002523"}],
    "references": [{"pubmedID": 32265444, "doi": "10.1038/s41467-020-15283-z"}],
}


def _mock_urlopen(response_data):
    """Helper to create a mock urlopen context manager."""
    mock_resp = MagicMock()
    mock_resp.read.return_value = json.dumps(response_data).encode()
    mock_resp.__enter__ = lambda s: s
    mock_resp.__exit__ = MagicMock(return_value=False)
    return mock_resp


@patch("qpx.core.http.urllib.request.urlopen")
def test_fetch_pride_metadata(mock_urlopen):
    """Basic fetch, no references, no organisms, 404 error, network error."""
    # Basic fetch
    mock_urlopen.return_value = _mock_urlopen(MOCK_PRIDE_RESPONSE)
    result = fetch_pride_metadata("PXD014414")
    assert result["project_title"] == "Quantitative proteomic landscape of metaplastic breast carcinoma"
    assert "proteomic analysis" in result["project_description"]
    assert result["pubmed_id"] == "32265444"
    assert result["doi"] == "10.6019/PXD014414"
    assert result["organisms"] == ["Homo sapiens (human)"]
    assert result["instruments"] == ["Q Exactive HF"]
    assert result["keywords"] == ["Proteomics", "Breast cancer"]

    # No references
    mock_urlopen.return_value = _mock_urlopen({**MOCK_PRIDE_RESPONSE, "references": []})
    assert fetch_pride_metadata("PXD014414")["pubmed_id"] is None

    # No organisms
    mock_urlopen.return_value = _mock_urlopen({**MOCK_PRIDE_RESPONSE, "organisms": None})
    assert fetch_pride_metadata("PXD014414")["organisms"] == []

    # 404 raises ValueError
    from urllib.error import HTTPError

    mock_urlopen.side_effect = HTTPError(url="", code=404, msg="Not Found", hdrs=None, fp=None)
    with pytest.raises(ValueError, match="not found"):
        fetch_pride_metadata("PXD999999")

    # Network error raises ConnectionError
    from urllib.error import URLError

    mock_urlopen.side_effect = URLError("Connection refused")
    with pytest.raises(ConnectionError, match="Cannot reach"):
        fetch_pride_metadata("PXD014414")


@patch("qpx.core.http.urllib.request.urlopen")
def test_enrich_from_pride(mock_urlopen, tmp_path):
    """enrich_from_pride updates dataset.parquet fields."""
    mock_urlopen.return_value = _mock_urlopen(MOCK_PRIDE_RESPONSE)

    from datetime import datetime

    from qpx.dataset import Dataset
    from qpx.writers.dataset import DatasetWriter

    dataset_path = tmp_path / "test.dataset.parquet"
    with DatasetWriter(dataset_path) as w:
        w.write_batch(
            [
                {
                    "project_accession": "PXD014414",
                    "creation_date": datetime.now().isoformat(),
                }
            ]
        )

    ds = Dataset(tmp_path)
    result = ds.enrich_from_pride()
    assert result["project_title"] is not None

    ds2 = Dataset(tmp_path)
    meta_df = ds2.dataset_meta.to_df()
    assert meta_df["project_title"].iloc[0] == "Quantitative proteomic landscape of metaplastic breast carcinoma"
    assert meta_df["pubmed_id"].iloc[0] == "32265444"
    ds.close()
    ds2.close()
