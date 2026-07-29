"""Plex-aware TMT channel labelling in the quantms feature adapter.

TMT10plex's 10th reporter is the standard ``131``; TMT11+ renames it ``131N`` to
pair with the added ``131C``. These lock that the feature adapter emits
``TMT131`` for a TMT10plex run (matching the CDAP converter and OpenMS/SDRF) and
``TMT131N``/``TMT131C`` for TMT11 -- so numeric-channel MSstats no longer mislabels
channel 10 as ``TMT131N`` on a TMT10 study.
"""

from __future__ import annotations

import pandas as pd

from qpx.converters.base import resolve_columns
from qpx.converters.quantms.feature_adapter import _FEATURE_MAP, QuantmsFeatureAdapter


def _adapter_with_channels(channels) -> QuantmsFeatureAdapter:
    """Build an adapter whose in-memory msstats table holds the given numeric channels."""
    adapter = QuantmsFeatureAdapter()
    values = ", ".join(f"({c})" for c in channels)
    adapter._conn.execute(f"CREATE TABLE msstats AS SELECT * FROM (VALUES {values}) AS t(Channel)")
    return adapter


def test_tmt10_channel_10_is_131():
    """A 10-channel run labels channel 10 as the standard TMT131, not TMT131N."""
    channel_map = _adapter_with_channels(range(1, 11))._resolve_tmt_channel_map()
    assert channel_map[1] == "TMT126"
    assert channel_map[10] == "TMT131"  # not TMT131N -- TMT10plex has no 131C


def test_tmt11_channel_10_is_131n_and_11_is_131c():
    """An 11-channel run keeps the TMT11 naming: channel 10 = TMT131N, 11 = TMT131C."""
    channel_map = _adapter_with_channels(range(1, 12))._resolve_tmt_channel_map()
    assert channel_map[10] == "TMT131N"
    assert channel_map[11] == "TMT131C"


def test_non_numeric_channels_leave_map_unchanged():
    """String channels (already labelled upstream) must not trigger the TMT10 relabel."""
    adapter = QuantmsFeatureAdapter()
    adapter._conn.execute("CREATE TABLE msstats AS SELECT * FROM (VALUES ('TMT131')) AS t(Channel)")
    assert adapter._resolve_tmt_channel_map()[10] == "TMT131N"


def test_tmt10_feature_labels_end_to_end():
    """A one-feature TMT10 batch yields TMT126..TMT131 with channel 10 = TMT131."""
    adapter = _adapter_with_channels(range(1, 11))
    adapter._tmt_channel_map = adapter._resolve_tmt_channel_map()
    adapter._modifications_meta = {}
    adapter._enzyme_name = None

    df = pd.DataFrame(
        {
            "PeptideSequence": ["PEPTIDE"] * 10,
            "ProteinName": ["P12345"] * 10,
            "Reference": ["run1.raw"] * 10,
            "Charge": [2] * 10,
            "Channel": list(range(1, 11)),
            "Intensity": [float(i) for i in range(1, 11)],
            "RetentionTime": [10.0] * 10,
        }
    )
    col_map = resolve_columns(_FEATURE_MAP, set(df.columns))
    records = adapter._transform_batch_isobaric(df, col_map, {}, {}, {}, "TMT")

    assert len(records) == 1
    labels = [entry["label"] for entry in records[0]["intensities"]]
    assert labels == [
        "TMT126",
        "TMT127N",
        "TMT127C",
        "TMT128N",
        "TMT128C",
        "TMT129N",
        "TMT129C",
        "TMT130N",
        "TMT130C",
        "TMT131",
    ]
