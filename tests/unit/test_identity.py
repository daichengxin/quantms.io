"""Tests for the derived mandatory identity ids (feature_id / psm_id / pg_id)."""

import pyarrow.parquet as pq

from qpx.core.data import FeatureSchema, PgSchema, PsmSchema
from qpx.core.data.identity import canonical, derive_id
from qpx.writers import FeatureWriter, PgWriter, PsmWriter
from tests.conftest import make_feature_record, make_pg_record, make_psm_record


def test_derive_id_is_deterministic():
    """Same input -> same id, on repeated calls."""
    values = ["PEPTIDEK", 2, "run_01", 120.5]
    assert derive_id(values) == derive_id(list(values))


def test_derive_id_distinguishes_inputs():
    """Different composites -> different ids (with overwhelming probability)."""
    assert derive_id(["PEPTIDEK", 2, "run_01", 120.5]) != derive_id(["PEPTIDEK", 3, "run_01", 120.5])
    assert derive_id(["PEPTIDEK", 2, "run_01", 120.5]) != derive_id(["PEPTIDEQ", 2, "run_01", 120.5])


def test_derive_id_handles_none():
    """None is encoded distinctly and does not collide with the string 'None'."""
    assert derive_id([None]) == derive_id([None])
    assert derive_id([None]) != derive_id(["None"])
    assert derive_id(["A", None, "B"]) != derive_id(["A", "B"])
    # int64 range
    for v in (derive_id([None]), derive_id(["A", None, "B"])):
        assert -(2**63) <= v < 2**63


def test_canonical_uses_separators():
    """The canonical encoding is unambiguous byte output."""
    assert canonical(["A", "B"]) == b"A\x1fB"
    assert canonical([None]) == b"\x00"


def test_schema_identity_metadata():
    """Each identity view exposes a single-column PK and its composite."""
    assert tuple(FeatureSchema._primary_key) == ("feature_id",)
    assert tuple(FeatureSchema.identity_composite) == ("peptidoform", "charge", "run_file_name", "rt")
    assert tuple(PsmSchema._primary_key) == ("psm_id",)
    assert tuple(PsmSchema.identity_composite) == ("peptidoform", "charge", "run_file_name", "scan")
    assert tuple(PgSchema._primary_key) == ("pg_id",)
    assert tuple(PgSchema.identity_composite) == ("anchor_protein", "grouped_runs", "label")


def test_feature_roundtrip_derives_id(tmp_path):
    """FeatureWriter stamps a non-null, unique feature_id matching derive_id."""
    path = tmp_path / "t.feature.parquet"
    records = [
        make_feature_record(sequence="PEPTIDEK", peptidoform="PEPTIDEK", charge=2, run_file_name="run_01"),
        make_feature_record(sequence="ANOTHERK", peptidoform="ANOTHERK", charge=3, run_file_name="run_01", scan=[2]),
    ]
    with FeatureWriter(path) as w:
        w.write_batch([dict(r) for r in records])

    table = pq.read_table(path)
    ids = table.column("feature_id").to_pylist()
    assert None not in ids
    assert len(set(ids)) == len(ids)
    for rec, got in zip(records, ids):
        expected = derive_id([rec["peptidoform"], rec["charge"], rec["run_file_name"], rec["rt"]])
        assert got == expected

    # Footer self-describes the key and composite.
    meta = {k.decode(): v.decode() for k, v in table.schema.metadata.items()}
    assert meta["primary_key"] == "feature_id"
    assert meta["identity_composite"] == "peptidoform,charge,run_file_name,rt"


def test_psm_roundtrip_derives_id(tmp_path):
    """PsmWriter stamps a non-null, unique psm_id matching derive_id."""
    path = tmp_path / "t.psm.parquet"
    records = [
        make_psm_record(sequence="PEPTIDEK", peptidoform="PEPTIDEK", charge=2, run_file_name="run_01", scan=[10]),
        make_psm_record(sequence="PEPTIDEK", peptidoform="PEPTIDEK", charge=2, run_file_name="run_01", scan=[11]),
    ]
    with PsmWriter(path) as w:
        w.write_batch([dict(r) for r in records])

    table = pq.read_table(path)
    ids = table.column("psm_id").to_pylist()
    assert None not in ids
    assert len(set(ids)) == len(ids)
    for rec, got in zip(records, ids):
        expected = derive_id([rec["peptidoform"], rec["charge"], rec["run_file_name"], rec["scan"]])
        assert got == expected


def test_pg_roundtrip_derives_id(tmp_path):
    """PgWriter stamps a non-null, unique pg_id matching derive_id of the flat composite."""
    path = tmp_path / "t.pg.parquet"
    records = [
        make_pg_record(anchor_protein="P1", run_file_name="run_01", intensities=[{"label": "L1", "intensity": 1.0}]),
        make_pg_record(anchor_protein="P2", run_file_name="run_01", intensities=[{"label": "L1", "intensity": 2.0}]),
    ]
    with PgWriter(path) as w:
        w.write_batch([dict(r) for r in records])

    table = pq.read_table(path)
    ids = table.column("pg_id").to_pylist()
    anchors = table.column("anchor_protein").to_pylist()
    grouped = table.column("grouped_runs").to_pylist()
    labels = table.column("label").to_pylist()
    assert None not in ids
    assert len(set(ids)) == len(ids)
    for i, got in enumerate(ids):
        assert got == derive_id([anchors[i], grouped[i], labels[i]])
