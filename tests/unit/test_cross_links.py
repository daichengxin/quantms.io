"""Cross-linking schema and round-trip tests."""

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from tests.conftest import make_psm_record


def test_cross_link_schema():
    """PSM schema includes cross_links as list<struct> with expected fields."""
    from qpx.core.data import PsmSchema

    field_names = [f.name for f in PsmSchema.get_arrow_schema()]
    assert "cross_links" in field_names

    field = PsmSchema.get_arrow_schema().field("cross_links")
    assert pa.types.is_list(field.type)
    assert pa.types.is_struct(field.type.value_type)

    xl_type = field.type.value_type
    xl_field_names = [xl_type.field(i).name for i in range(xl_type.num_fields)]
    for expected in (
        "xl_type",
        "partner_sequence",
        "partner_peptidoform",
        "donor_position",
        "acceptor_position",
        "linker_name",
        "linker_accession",
        "linker_mass",
    ):
        assert expected in xl_field_names


def test_cross_link_round_trip(tmp_path):
    """Round-trip: inter XL, null XL, dead-end XL, multiple XLs."""
    from qpx.writers import PsmWriter

    # Inter cross-link
    record = make_psm_record()
    record["cross_links"] = [
        {
            "xl_type": "inter",
            "partner_sequence": "PEPTIDEK",
            "partner_peptidoform": "PEPTIDEK",
            "donor_position": 5,
            "acceptor_position": 3,
            "linker_name": "DSS",
            "linker_accession": "XLMOD:02001",
            "linker_mass": 138.0681,
        }
    ]
    out = tmp_path / "xl.psm.parquet"
    with PsmWriter(out) as w:
        w.write_batch([record])
    table = pq.read_table(out)
    xl = table.column("cross_links")[0].as_py()
    assert len(xl) == 1
    assert xl[0]["xl_type"] == "inter"
    assert xl[0]["linker_mass"] == pytest.approx(138.0681)

    # Null cross_links for regular PSM
    record2 = make_psm_record()
    out2 = tmp_path / "regular.psm.parquet"
    with PsmWriter(out2) as w:
        w.write_batch([record2])
    assert pq.read_table(out2).column("cross_links")[0].as_py() is None

    # Dead-end cross-link
    record3 = make_psm_record()
    record3["cross_links"] = [
        {
            "xl_type": "dead-end",
            "partner_sequence": None,
            "partner_peptidoform": None,
            "donor_position": 5,
            "acceptor_position": None,
            "linker_name": "DSS",
            "linker_accession": None,
            "linker_mass": 156.0786,
        }
    ]
    out3 = tmp_path / "de.psm.parquet"
    with PsmWriter(out3) as w:
        w.write_batch([record3])
    xl = pq.read_table(out3).column("cross_links")[0].as_py()
    assert xl[0]["xl_type"] == "dead-end"
    assert xl[0]["partner_sequence"] is None

    # Multiple cross-links
    record4 = make_psm_record()
    record4["cross_links"] = [
        {
            "xl_type": "inter",
            "partner_sequence": "PEPTIDEK",
            "partner_peptidoform": "PEPTIDEK",
            "donor_position": 5,
            "acceptor_position": 3,
            "linker_name": "DSSO",
            "linker_accession": None,
            "linker_mass": 158.0038,
        },
        {
            "xl_type": "intra",
            "partner_sequence": None,
            "partner_peptidoform": None,
            "donor_position": 2,
            "acceptor_position": 7,
            "linker_name": "DSSO",
            "linker_accession": None,
            "linker_mass": 158.0038,
        },
    ]
    out4 = tmp_path / "multi_xl.psm.parquet"
    with PsmWriter(out4) as w:
        w.write_batch([record4])
    xl = pq.read_table(out4).column("cross_links")[0].as_py()
    assert len(xl) == 2
    assert xl[0]["xl_type"] == "inter"
    assert xl[1]["xl_type"] == "intra"
