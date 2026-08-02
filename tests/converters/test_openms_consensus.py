"""consensusXML -> QPX converter (interim path).

Reads a tiny literal consensusXML fixture (no pyopenms construction) so the
extraction runs in CI without a large fixture file and without depending on the
pyopenms setter APIs, which differ across versions.
"""

import duckdb
import pytest

pytest.importorskip("pyopenms")
import pyopenms as oms  # noqa: E402

from qpx.converters.openms_consensus.converter import OpenMSConsensusConverter  # noqa: E402
from qpx.converters.openms_consensus.feature_adapter import to_proforma  # noqa: E402


@pytest.mark.parametrize(
    ("openms_seq", "expected"),
    [
        ("PEPTIDEK", "PEPTIDEK"),
        (".(TMT6plex)THSQEEM(Oxidation)QHMQR", "[UNIMOD:737]-THSQEEM[UNIMOD:35]QHMQR"),
        ("C(Carbamidomethyl)PEPTIDEK", "C[UNIMOD:4]PEPTIDEK"),
        ("PEPTIDER.(Amidated)", "PEPTIDER-[UNIMOD:2]"),
    ],
)
def test_to_proforma(openms_seq, expected):
    assert to_proforma(oms.AASequence.fromString(openms_seq)) == expected


# A 2-channel isobaric consensusXML written as a literal fixture: one peptide,
# 2 TMT channels (126/127), 1 protein. Kept as text — not built through pyopenms
# constructors — so the test exercises only the *read* path and is immune to
# pyopenms-version drift in the setter APIs (e.g. list vs PeptideIdentificationList).
_TMT_CONSENSUSXML = """<?xml version="1.0" encoding="ISO-8859-1"?>
<consensusXML version="1.7" experiment_type="labeled_MS2" xsi:noNamespaceSchemaLocation="https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/SCHEMAS/ConsensusXML_1_7.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
	<IdentificationRun id="PI_0" date="0000-00-00T00:00:00" search_engine="" search_engine_version="">
		<SearchParameters db="" db_version="" taxonomy="" mass_type="monoisotopic" charges="" enzyme="unknown_enzyme" missed_cleavages="0" precursor_peak_tolerance="0" precursor_peak_tolerance_ppm="false" peak_mass_tolerance="0" peak_mass_tolerance_ppm="false" >
		</SearchParameters>
		<ProteinIdentification score_type="" higher_score_better="true" significance_threshold="0">
			<ProteinHit id="PH_0" accession="P12345" score="0" sequence="">
			</ProteinHit>
		</ProteinIdentification>
	</IdentificationRun>
	<mapList count="2">
		<map id="0" name="run_01.mzML" unique_id="1" label="tmt6plex_126" size="1">
		</map>
		<map id="1" name="run_01.mzML" unique_id="2" label="tmt6plex_127" size="1">
		</map>
	</mapList>
	<consensusElementList>
		<consensusElement id="e_0" quality="0.0" charge="2">
			<centroid rt="100.0" mz="450.25" it="0.0"/>
			<groupedElementList>
				<element map="0" id="0" rt="100.0" mz="450.25" it="1000.0"/>
				<element map="1" id="1" rt="100.0" mz="450.25" it="2000.0"/>
			</groupedElementList>
			<PeptideIdentification identification_run_ref="PI_0" score_type="" higher_score_better="true" significance_threshold="0" MZ="450.26" RT="100" spectrum_reference="controllerType=0 controllerNumber=1 scan=42" >
				<PeptideHit score="0" sequence="PEPTIDEK" charge="2" protein_refs="PH_0">
					<UserParam type="string" name="target_decoy" value="target"/>
					<UserParam type="float" name="Posterior Error Probability_score" value="1.0e-03"/>
				</PeptideHit>
			</PeptideIdentification>
		</consensusElement>
	</consensusElementList>
</consensusXML>
"""


def _write_tmt_consensusxml(path):
    """Write the literal 2-channel TMT consensusXML fixture to ``path``."""
    path.write_text(_TMT_CONSENSUSXML)


def test_converter_rejects_invalid_structure_before_creating_output(tmp_path):
    out = tmp_path / "out"
    with pytest.raises(ValueError, match="Unknown.*bogus"):
        OpenMSConsensusConverter().convert(
            str(tmp_path / "missing.consensusXML"),
            str(out),
            structures=("bogus",),
        )
    assert not out.exists()


def test_converter_requires_sdrf_for_requested_metadata(tmp_path):
    out = tmp_path / "out"
    with pytest.raises(ValueError, match="SDRF.*run"):
        OpenMSConsensusConverter().convert(
            str(tmp_path / "missing.consensusXML"),
            str(out),
            structures=("run",),
        )
    assert not out.exists()


def test_converter_writes_only_requested_sdrf_structure(tmp_path):
    sdrf = tmp_path / "test.sdrf.tsv"
    sdrf.write_text(
        "source name\tcharacteristics[organism]\tcharacteristics[organism part]\t"
        "comment[data file]\tcomment[label]\n"
        "S1\tHomo sapiens\tliver\trun_01.raw\tlabel free sample\n"
    )
    out = tmp_path / "out"
    written = OpenMSConsensusConverter().convert(
        str(tmp_path / "unused.consensusXML"),
        str(out),
        sdrf_path=str(sdrf),
        structures=("run",),
    )

    assert set(written) == {"run"}
    assert written["run"].exists()
    assert not (out / "openms.sample.parquet").exists()


def test_pg_peptide_counts_are_per_protein(monkeypatch):
    from qpx.converters.openms_consensus import pg_adapter

    class Header:
        filename = "run_01.mzML"
        label = ""

    class ConsensusMap:
        @staticmethod
        def getColumnHeaders():
            return {0: Header()}

        @staticmethod
        def getExperimentType():
            return "label-free"

        @staticmethod
        def getProteinIdentifications():
            return [object()]

    monkeypatch.setattr(
        pg_adapter,
        "_protein_maps",
        lambda _cm: (
            {"P1": {"PEPA", "PEPB"}, "P2": {"PEPB"}},
            {"P1": {"run_01"}, "P2": {"run_01"}},
            {"P1": {("PEPA", 2), ("PEPB", 2)}, "P2": {("PEPB", 2)}},
        ),
    )
    monkeypatch.setattr(pg_adapter, "_protein_hit_meta", lambda _prot: ({}, {}, {}))
    monkeypatch.setattr(pg_adapter, "_build_groups", lambda _prot: [["P1", "P2"]])

    records = pg_adapter.consensus_protein_groups_to_records(cm=ConsensusMap())

    assert records[0]["peptide_counts"]["unique_sequences"] == 2
    assert records[0]["peptides"] == [
        {"protein_name": "P1", "peptide_count": 2},
        {"protein_name": "P2", "peptide_count": 1},
    ]


def test_pg_uses_id_merge_index_after_unmapped_map_index():
    from qpx.converters.openms_consensus.pg_adapter import _protein_maps

    class Header:
        def __init__(self, filename):
            self.filename = filename

    class ConsensusMap:
        @staticmethod
        def getColumnHeaders():
            return {0: Header("run_01.mzML"), 1: Header("run_02.mzML")}

        @staticmethod
        def getUnassignedPeptideIdentifications():
            return [PeptideIdentification()]

        def __iter__(self):
            return iter(())

    class Sequence:
        @staticmethod
        def toUnmodifiedString():
            return "PEPTIDE"

        @staticmethod
        def toUniModString():
            return "PEPTIDE"

    class Evidence:
        @staticmethod
        def getProteinAccession():
            return "P1"

    class Hit:
        @staticmethod
        def getSequence():
            return Sequence()

        @staticmethod
        def getCharge():
            return 2

        @staticmethod
        def getPeptideEvidences():
            return [Evidence()]

    class PeptideIdentification:
        values = {"map_index": 99, "id_merge_index": 1}

        def metaValueExists(self, key):
            return key in self.values

        def getMetaValue(self, key):
            return self.values[key]

        @staticmethod
        def getHits():
            return [Hit()]

    peptides, runs, features = _protein_maps(ConsensusMap())
    assert peptides == {"P1": {"PEPTIDE"}}
    assert runs == {"P1": {"run_02"}}
    assert features == {"P1": {("PEPTIDE", 2)}}


def test_consensusxml_to_qpx_feature_has_channels_pg_is_identification_only(tmp_path):
    cx = tmp_path / "test.consensusXML"
    _write_tmt_consensusxml(cx)
    out = tmp_path / "out"
    written = OpenMSConsensusConverter().convert(str(cx), str(out), output_prefix="t", structures=("feature", "psm", "pg"))
    con = duckdb.connect()

    feat = con.execute(
        f"SELECT peptidoform, charge, run_file_name, intensities FROM read_parquet('{written['feature']}')"
    ).fetchall()
    assert len(feat) == 1
    pep, charge, run, intensities = feat[0]
    assert pep == "PEPTIDEK" and charge == 2 and run == "run_01"
    labels = {e["label"]: e["intensity"] for e in intensities}
    assert labels == {"TMT126": 1000.0, "TMT127": 2000.0}  # both channels, canonicalized, quant kept

    psm = con.execute(f"SELECT peptidoform, scan FROM read_parquet('{written['psm']}')").fetchall()
    assert psm and psm[0][0] == "PEPTIDEK" and list(psm[0][1]) == [42]

    pg = con.execute(f"SELECT anchor_protein, label, intensity FROM read_parquet('{written['pg']}')").fetchall()
    assert pg and all(anchor == "P12345" for anchor, _, _ in pg)
    # interim: one row per channel (the quantification slot), intensity null until
    # OpenMS -out_qpx fills it. The 2 TMT channels each get a placeholder row.
    assert {label for _, label, _ in pg} == {"TMT126", "TMT127"}
    assert all(intensity is None for _, _, intensity in pg)


def test_label_free_consensusxml_uses_lfq_labels(tmp_path):
    """OpenMS spells its label-free experiment type with a hyphen."""
    cx = tmp_path / "label_free.consensusXML"
    cx.write_text(_TMT_CONSENSUSXML.replace('experiment_type="labeled_MS2"', 'experiment_type="label-free"'))
    out = tmp_path / "out"
    written = OpenMSConsensusConverter().convert(
        str(cx),
        str(out),
        output_prefix="lfq",
        structures=("feature", "pg"),
    )
    con = duckdb.connect()
    feature_labels = con.execute(
        "SELECT UNNEST(intensities).label FROM read_parquet($1)",
        [str(written["feature"])],
    ).fetchall()
    pg_labels = con.execute(
        "SELECT label FROM read_parquet($1)",
        [str(written["pg"])],
    ).fetchall()
    con.close()

    assert {label for (label,) in feature_labels} == {"LFQ"}
    assert {label for (label,) in pg_labels} == {"LFQ"}
