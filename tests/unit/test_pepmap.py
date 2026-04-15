"""Peptide-protein mapping tests."""

import pyarrow.parquet as pq


def make_mapping_record(**overrides):
    defaults = {
        "sequence": "PEPTIDEK",
        "peptidoform": "PEPTIDEK",
        "pg_accessions": [
            {"accession": "P12345", "start": 100, "end": 107, "pre": "K", "post": "R"},
        ],
        "is_unique": True,
    }
    defaults.update(overrides)
    return defaults


def test_pepmap_schema_and_round_trip(tmp_path):
    """Schema loads; round-trip writes and reads; nullable fields work."""
    from qpx.core.data import PepMapSchema
    from qpx.writers import PepMapWriter

    # Schema loads with expected fields
    field_names = [f.name for f in PepMapSchema.get_arrow_schema()]
    assert "sequence" in field_names
    assert "peptidoform" in field_names
    assert "pg_accessions" in field_names
    assert "is_unique" in field_names

    # Round trip
    records = [
        make_mapping_record(),
        make_mapping_record(
            sequence="ANOTHERK",
            peptidoform="ANOTHERK",
            pg_accessions=[
                {
                    "accession": "P67890",
                    "start": None,
                    "end": None,
                    "pre": None,
                    "post": None,
                },
            ],
            is_unique=False,
        ),
    ]
    out = tmp_path / "test.pepmap.parquet"
    with PepMapWriter(out) as w:
        w.write_batch(records)
    table = pq.read_table(out)
    assert table.num_rows == 2
    assert table.column("sequence")[0].as_py() == "PEPTIDEK"
    assert table.column("peptidoform")[1].as_py() == "ANOTHERK"

    # Nullable fields
    record = make_mapping_record(pg_accessions=None, is_unique=None)
    out2 = tmp_path / "nullable.pepmap.parquet"
    with PepMapWriter(out2) as w:
        w.write_batch([record])
    assert pq.read_table(out2).num_rows == 1


def test_pepmap_dataset_integration(tmp_path, dataset_dir):
    """Dataset discovers pepmap; filter_by_protein; save_structure round-trip."""
    import qpx
    from qpx.writers import PepMapWriter

    # Create pepmap with multiple records
    records = [
        make_mapping_record(
            sequence="PEP1K",
            peptidoform="PEP1K",
            pg_accessions=[
                {
                    "accession": "P12345",
                    "start": None,
                    "end": None,
                    "pre": None,
                    "post": None,
                }
            ],
        ),
        make_mapping_record(
            sequence="PEP2K",
            peptidoform="PEP2K",
            pg_accessions=[
                {
                    "accession": "P67890",
                    "start": None,
                    "end": None,
                    "pre": None,
                    "post": None,
                }
            ],
        ),
        make_mapping_record(
            sequence="PEP3K",
            peptidoform="PEP3K",
            pg_accessions=[
                {
                    "accession": "P12345",
                    "start": None,
                    "end": None,
                    "pre": None,
                    "post": None,
                },
                {
                    "accession": "P67890",
                    "start": None,
                    "end": None,
                    "pre": None,
                    "post": None,
                },
            ],
            is_unique=False,
        ),
    ]
    out = tmp_path / "exp.pepmap.parquet"
    with PepMapWriter(out) as w:
        w.write_batch(records)

    # Dataset discovers pepmap
    ds = qpx.open_dataset(str(tmp_path))
    assert ds.pepmap is not None
    assert ds.pepmap.count() == 3

    # Filter by protein
    filtered = ds.pepmap.by_protein("P12345")
    assert filtered.count() == 2
    ds.close()

    # save_structure round-trip
    ds3 = qpx.open_dataset(str(dataset_dir))
    save_records = [
        make_mapping_record(),
        make_mapping_record(sequence="ANOTHERK", peptidoform="ANOTHERK"),
    ]
    path = ds3.save_structure(save_records, "pepmap")
    assert path.exists()
    ds3.refresh()
    assert ds3.pepmap is not None
    assert ds3.pepmap.count() == 2
    ds3.close()
