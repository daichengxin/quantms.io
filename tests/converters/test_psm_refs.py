"""Tests for reading authoritative PSM reference fields from a consensusXML."""

from __future__ import annotations

import textwrap

from qpx.converters.channel_labels import parse_consensusxml_maplist
from qpx.converters.openms.psm_refs import psm_references_from_consensusxml

# Minimal consensusXML: two maps (rep1/rep2), one protein, one assigned PSM on
# map 1 (rep2) and one unassigned decoy PSM on map 0 (rep1) that share scan=100.
_CXML = textwrap.dedent(
    """\
    <?xml version="1.0" encoding="ISO-8859-1"?>
    <consensusXML version="1.7">
      <mapList count="2">
        <map id="0" name="/data/pool_rep1.mzML" label="" size="10"/>
        <map id="1" name="/data/pool_rep2.mzML" label="" size="10"/>
      </mapList>
      <consensusElementList>
        <consensusElement id="e_0">
          <centroid rt="600.0" mz="500.25" intensity="1000"/>
          <PeptideIdentification identification_run_ref="PI_0" score_type="q-value"
              higher_score_better="false" MZ="500.25" RT="600.0"
              spectrum_reference="controllerType=0 controllerNumber=1 scan=100" >
            <PeptideHit score="0.01" sequence="PEPT(Phospho)IDEK" charge="2" protein_refs="PH_0">
              <UserParam type="string" name="target_decoy" value="target"/>
              <UserParam type="float" name="percolator_score" value="-1.234567890123e-03"/>
              <UserParam type="float" name="Posterior Error Probability_score" value="0.2"/>
            </PeptideHit>
            <UserParam type="int" name="map_index" value="1"/>
            <UserParam type="string" name="feature_id" value="9999888877776666"/>
          </PeptideIdentification>
        </consensusElement>
      </consensusElementList>
      <UnassignedPeptideIdentification identification_run_ref="PI_0" score_type="q-value"
          higher_score_better="false" MZ="640.30" RT="720.0"
          spectrum_reference="controllerType=0 controllerNumber=1 scan=100" >
        <PeptideHit score="0.05" sequence="DECOYPEPK" charge="3" protein_refs="PH_0">
          <UserParam type="string" name="target_decoy" value="decoy"/>
          <UserParam type="float" name="percolator_score" value="-7.6e-04"/>
          <UserParam type="float" name="Posterior Error Probability_score" value="0.6"/>
        </PeptideHit>
        <UserParam type="int" name="map_index" value="0"/>
        <UserParam type="string" name="feature_id" value="not mapped"/>
      </UnassignedPeptideIdentification>
    </consensusXML>
    """
)


def _write(tmp_path):
    # ProteinHit lives in the ProteinIdentification section; inject one so
    # protein_refs resolve to an accession.
    cxml = _CXML.replace(
        "<consensusElementList>",
        '<ProteinIdentification><ProteinHit id="PH_0" accession="sp|P12345|TEST"/>'
        "</ProteinIdentification>\n  <consensusElementList>",
    )
    path = tmp_path / "test.consensusXML"
    path.write_text(cxml)
    return str(path)


def test_reads_run_file_from_map_index(tmp_path):
    path = _write(tmp_path)
    refs = psm_references_from_consensusxml(path, parse_consensusxml_maplist(path))
    assert len(refs) == 2
    by_seq = {r.sequence: r for r in refs}
    assigned = by_seq["PEPT(Phospho)IDEK"]
    unassigned = by_seq["DECOYPEPK"]
    # map_index=1 -> rep2 (NOT the first run), which is exactly the bug we fix.
    assert assigned.run_file_name == "/data/pool_rep2.mzML"
    assert unassigned.run_file_name == "/data/pool_rep1.mzML"


def test_extracts_scan_scores_decoy_proteins_and_feature_link(tmp_path):
    path = _write(tmp_path)
    refs = psm_references_from_consensusxml(path, parse_consensusxml_maplist(path))
    assigned = next(r for r in refs if r.sequence == "PEPT(Phospho)IDEK")
    unassigned = next(r for r in refs if r.sequence == "DECOYPEPK")

    assert assigned.scan == 100 and unassigned.scan == 100
    assert assigned.charge == 2 and assigned.observed_mz == 500.25
    assert assigned.percolator_score == -1.234567890123e-03
    assert assigned.posterior_error_probability == 0.2
    assert assigned.is_decoy is False and unassigned.is_decoy is True
    assert assigned.protein_accessions == ["sp|P12345|TEST"]
    # assigned carries the OpenMS feature id; unassigned's "not mapped" -> None
    assert assigned.assigned is True and assigned.openms_feature_id == "9999888877776666"
    assert unassigned.assigned is False and unassigned.openms_feature_id is None


def test_missing_file_returns_empty(tmp_path):
    assert psm_references_from_consensusxml(str(tmp_path / "nope.consensusXML"), {}) == []
