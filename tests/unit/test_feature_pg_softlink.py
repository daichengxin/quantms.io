"""Unit tests for the computed feature<->pg softlink (bigbio/qpx#269).

The feature->pg association is derived on read by ``Dataset.link_feature_pg()``,
a label-aware join of the ``feature`` and ``pg`` views. These tests build tiny
synthetic datasets directly so the label-awareness and the identified-but-not-
quantified behaviour can be asserted in isolation, without any converter.
"""

from __future__ import annotations

from qpx.dataset import Dataset
from qpx.writers import FeatureWriter, PgWriter
from tests.conftest import make_feature_record, make_pg_record


def _feature(sequence, run, labels, membership, rt):
    """A feature carrying ``labels`` (list of channel names) for ``membership``."""
    rec = make_feature_record(
        sequence=sequence,
        peptidoform=sequence,
        run_file_name=run,
        intensities=[{"label": lbl, "intensity": 10.0 * (i + 1)} for i, lbl in enumerate(labels)],
    )
    rec["rt"] = rt
    rec["pg_accessions"] = [{"accession": acc, "start": None, "end": None} for acc in membership]
    rec["anchor_protein"] = membership[0]
    return rec


def _write(ds_dir, features, pgs):
    ds_dir.mkdir(parents=True, exist_ok=True)
    with FeatureWriter(ds_dir / "exp.feature.parquet") as w:
        w.write_batch(features)
    with PgWriter(ds_dir / "exp.pg.parquet") as w:
        w.write_batch(pgs)


def test_softlink_is_label_aware_no_over_link(tmp_path):
    """#269: a same-membership pg row that exists only for a label the feature
    lacks must NOT be linked. Feature carries channel 'L'; pg has both 'L' and 'H'
    for the same membership + run -> the feature links to 'L' only."""
    ds_dir = tmp_path / "labelaware"
    feature = _feature("PEPTIDEF", "run_01", ["L"], ["P1", "P2"], rt=100.0)
    # One pg record with two labels -> the writer explodes it into two rows
    # (label L and label H), same membership + grouped_runs.
    pg = make_pg_record(
        pg_accessions=["P1", "P2"],
        run_file_name="run_01",
        intensities=[{"label": "L", "intensity": 100.0}, {"label": "H", "intensity": 50.0}],
    )
    _write(ds_dir, [feature], [pg])

    with Dataset(ds_dir) as ds:
        link = ds.link_feature_pg().fetchall()
        pg_rows = ds.sql("SELECT pg_id, label FROM pg").fetchall()

    id_by_label = {label: pg_id for pg_id, label in pg_rows}
    assert set(id_by_label) == {"L", "H"}

    linked_labels = {label for _fid, _pid, label in link}
    linked_ids = {pid for _fid, pid, _label in link}
    assert linked_labels == {"L"}, "feature must link only to the label it carries"
    assert id_by_label["L"] in linked_ids
    assert id_by_label["H"] not in linked_ids, "must not over-link to the missing channel (#269)"


def test_softlink_identified_but_not_quantified_has_no_link(tmp_path):
    """A feature whose group has no matching pg row produces no softlink edge
    (identified but not quantified) — it is not an error."""
    ds_dir = tmp_path / "idonly"
    quantified = _feature("PEPTIDEA", "run_01", ["L"], ["P1", "P2"], rt=100.0)
    id_only = _feature("PEPTIDEZ", "run_01", ["L"], ["P9"], rt=200.0)  # P9 has no pg row
    pg = make_pg_record(
        pg_accessions=["P1", "P2"],
        run_file_name="run_01",
        intensities=[{"label": "L", "intensity": 100.0}],
    )
    _write(ds_dir, [quantified, id_only], [pg])

    with Dataset(ds_dir) as ds:
        link = ds.link_feature_pg().fetchall()
        unlinked = {row[0] for row in ds.features_without_pg_link().fetchall()}
        feat_ids = {row[0]: row[1] for row in ds.sql("SELECT feature_id, sequence FROM feature").fetchall()}

    seq_by_id = feat_ids
    linked_ids = {fid for fid, _pid, _label in link}
    quantified_id = next(fid for fid, seq in seq_by_id.items() if seq == "PEPTIDEA")
    idonly_id = next(fid for fid, seq in seq_by_id.items() if seq == "PEPTIDEZ")

    assert quantified_id in linked_ids
    assert idonly_id not in linked_ids
    assert idonly_id in unlinked
    assert quantified_id not in unlinked
