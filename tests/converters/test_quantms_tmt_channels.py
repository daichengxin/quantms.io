"""Plex-aware TMT channel labelling in the quantms feature adapter.

TMT10plex's 10th reporter is the standard ``131``; TMT11+ renames it ``131N`` to
pair with the added ``131C``. These lock that the feature adapter emits
``TMT131`` for a TMT10plex run (matching OpenMS/SDRF and the shared channel_map)
and ``TMT131N``/``TMT131C`` for TMT11 -- so numeric-channel MSstats no longer
mislabels channel 10 as ``TMT131N`` on a TMT10 study.

Channel labels are resolved through the shared vocabulary
(:mod:`qpx.converters.channel_labels`); these tests exercise it through the
feature adapter's isobaric transform.
"""

from __future__ import annotations

import pandas as pd

from qpx.converters.base import resolve_columns
from qpx.converters.channel_labels import resolve_channel_labels
from qpx.converters.quantms.feature_adapter import _FEATURE_MAP, QuantmsFeatureAdapter


def _transform(channels, experiment_type="TMT"):
    """Run the isobaric transform on a one-feature batch with the given channels."""
    adapter = QuantmsFeatureAdapter()
    adapter._modifications_meta = {}
    adapter._enzyme_name = None
    n = len(channels)
    df = pd.DataFrame(
        {
            "PeptideSequence": ["PEPTIDE"] * n,
            "ProteinName": ["P12345"] * n,
            "Reference": ["run1.raw"] * n,
            "Charge": [2] * n,
            "Channel": list(channels),
            "Intensity": [float(i + 1) for i in range(n)],
            "RetentionTime": [10.0] * n,
        }
    )
    col_map = resolve_columns(_FEATURE_MAP, set(df.columns))
    records = adapter._transform_batch_isobaric(df, col_map, {}, {}, {}, experiment_type)
    assert len(records) == 1
    return [e["label"] for e in records[0]["intensities"]]


def test_tmt10_channel_10_is_131():
    """A 10-channel run labels channel 10 as the standard TMT131, not TMT131N."""
    labels = resolve_channel_labels("TMT", None, range(1, 11))
    assert labels[1] == "TMT126"
    assert labels[10] == "TMT131"  # not TMT131N -- TMT10plex has no 131C


def test_tmt11_channel_10_is_131n_and_11_is_131c():
    """An 11-channel run keeps the TMT11 naming: channel 10 = TMT131N, 11 = TMT131C."""
    labels = resolve_channel_labels("TMT", None, range(1, 12))
    assert labels[10] == "TMT131N"
    assert labels[11] == "TMT131C"


def test_tmt10_feature_labels_end_to_end():
    """A one-feature TMT10 batch yields TMT126..TMT131 with channel 10 = TMT131."""
    assert _transform(range(1, 11)) == [
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


def test_tmt11_feature_labels_end_to_end():
    """A one-feature TMT11 batch yields channel 10 = TMT131N, 11 = TMT131C."""
    labels = _transform(range(1, 12))
    assert labels[9] == "TMT131N"
    assert labels[10] == "TMT131C"


def test_non_numeric_channel_label_passes_through():
    """A channel that is already a reporter string is preserved, not re-indexed."""
    assert _transform(["TMT126"]) == ["TMT126"]


def test_mixed_tmt10_tmt11_batch_labels_each_run_by_its_own_plex():
    """A batch combining a TMT10 run and a TMT11 run labels each by its own plex.

    Channel labels are resolved per run, so the TMT10 run's channel 10 is TMT131
    while the TMT11 run's channel 10 is TMT131N -- a single batch-wide resolution
    would see the union (1..11) and wrongly label both channel 10 as TMT131N.
    """
    adapter = QuantmsFeatureAdapter()
    adapter._modifications_meta = {}
    adapter._enzyme_name = None
    # No SDRF labels -> plex is inferred from each run's channel count.
    rows = []
    for ch in range(1, 11):  # run_a: TMT10
        rows.append(("PEPA", "P1", "run_a.raw", 2, ch, float(ch), 10.0))
    for ch in range(1, 12):  # run_b: TMT11
        rows.append(("PEPB", "P2", "run_b.raw", 2, ch, float(ch), 10.0))
    df = pd.DataFrame(
        rows,
        columns=["PeptideSequence", "ProteinName", "Reference", "Charge", "Channel", "Intensity", "RetentionTime"],
    )
    col_map = resolve_columns(_FEATURE_MAP, set(df.columns))
    records = adapter._transform_batch_isobaric(df, col_map, {}, {}, {}, "TMT")
    by_run = {r["run_file_name"]: [e["label"] for e in r["intensities"]] for r in records}
    assert by_run["run_a"][9] == "TMT131"  # TMT10 channel 10
    assert "TMT131N" not in by_run["run_a"]
    assert by_run["run_b"][9] == "TMT131N"  # TMT11 channel 10
    assert by_run["run_b"][10] == "TMT131C"


def test_sdrf_itraq_family_overrides_numeric_channel_fallback():
    """Numeric 1..4 channels are iTRAQ when the SDRF declares iTRAQ labels."""
    with QuantmsFeatureAdapter() as adapter:
        adapter._conn.execute('CREATE TABLE msstats ("Channel" INTEGER)')
        adapter._conn.execute("INSERT INTO msstats VALUES (1), (2), (3), (4)")
        assert adapter._resolve_experiment_type({"ITRAQ114", "ITRAQ115", "ITRAQ116", "ITRAQ117"}) == "iTRAQ"
