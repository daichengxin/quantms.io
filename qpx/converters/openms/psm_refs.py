"""Authoritative PSM reference fields from an OpenMS consensusXML.

OpenMS ``ProteomicsLFQ -out_qpx`` writes a ``psm.parquet`` whose per-PSM
``run_file_name`` is wrong (every row is stamped with the first run's file —
the ``map_index`` is dropped; OpenMS#9872) and which duplicates a handful of
spectra (one feature-attached + one unassigned copy; OpenMS#9871). The
companion consensusXML carries the correct information, so this module reads the
reference fields from there and the converter overlays them onto the parquet's
(authoritative) chemistry, matched per PSM on the exact key
``(plain_sequence, charge, scan)`` — all bit-stable values, so the pairing is
1:1 without relying on document order or drift-prone floats.

Each ``PeptideHit`` is emitted as one :class:`PsmRef`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from defusedxml import DefusedXmlException
from defusedxml.ElementTree import ParseError, iterparse

# The consensusXML score written by quantms' Percolator step; full-precision, so
# it discriminates PSM rows between the consensusXML and the parquet.
_MATCH_SCORE_PARAM = "percolator_score"
_PEP_PARAM = "Posterior Error Probability_score"


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _scan_of(spectrum_reference: Optional[str]) -> Optional[int]:
    """Extract the integer scan from an OpenMS ``spectrum_reference`` string."""
    for token in (spectrum_reference or "").split():
        if token.startswith("scan="):
            try:
                return int(token.split("=", 1)[1])
            except ValueError:
                return None
    return None


@dataclass
class PsmRef:  # pylint: disable=too-many-instance-attributes
    """One consensusXML PeptideHit's reference fields (not its chemistry)."""

    run_file_name: Optional[str]
    scan: Optional[int]
    charge: Optional[int]
    percolator_score: Optional[float]
    posterior_error_probability: Optional[float]
    observed_mz: Optional[float]
    rt: Optional[float]
    is_decoy: bool
    protein_accessions: list[str] = field(default_factory=list)
    openms_feature_id: Optional[str] = None
    assigned: bool = False
    sequence: Optional[str] = None


def _user_params(element) -> dict[str, str]:
    return {child.attrib.get("name", ""): child.attrib.get("value", "") for child in element if _local(child.tag) == "UserParam"}


def _hit_ref(pi_element, maplist: dict[int, dict[str, str]], proteins: dict[str, str]) -> Optional[PsmRef]:
    """Build a :class:`PsmRef` from a (Unassigned)PeptideIdentification element."""
    hits = [c for c in pi_element if _local(c.tag) == "PeptideHit"]
    if not hits:
        return None
    hit = hits[0]  # rank-1; OpenMS -out_qpx writes one hit per ID
    id_params = _user_params(pi_element)
    hit_params = _user_params(hit)

    map_index = id_params.get("map_index")
    run_file_name = None
    if map_index is not None:
        try:
            run_file_name = (maplist.get(int(map_index)) or {}).get("name") or None
        except ValueError:
            run_file_name = None

    def _float(value: Optional[str]) -> Optional[float]:
        try:
            return float(value) if value is not None else None
        except ValueError:
            return None

    accessions = [proteins.get(ref, ref) for ref in (hit.attrib.get("protein_refs") or "").split()]
    feature_id = id_params.get("feature_id")
    if feature_id in (None, "", "not mapped"):
        feature_id = None

    return PsmRef(
        run_file_name=run_file_name,
        scan=_scan_of(pi_element.attrib.get("spectrum_reference")),
        charge=int(hit.attrib["charge"]) if hit.attrib.get("charge") else None,
        percolator_score=_float(hit_params.get(_MATCH_SCORE_PARAM)),
        posterior_error_probability=_float(hit_params.get(_PEP_PARAM)),
        observed_mz=_float(pi_element.attrib.get("MZ")),
        rt=_float(pi_element.attrib.get("RT")),
        is_decoy=hit_params.get("target_decoy") == "decoy",
        protein_accessions=accessions,
        openms_feature_id=feature_id,
        assigned=_local(pi_element.tag) == "PeptideIdentification",
        sequence=hit.attrib.get("sequence"),
    )


def psm_references_from_consensusxml(
    consensusxml_path: str,
    maplist: dict[int, dict[str, str]],
) -> list[PsmRef]:
    """Read every PSM's reference fields from a consensusXML.

    Emits one :class:`PsmRef` per ``PeptideHit`` across both assigned
    (``PeptideIdentification``) and unassigned (``UnassignedPeptideIdentification``)
    identifications. ``maplist`` is the shared ``parse_consensusxml_maplist``
    result (``map_index`` -> ``{"name": file, ...}``). Returns ``[]`` on a missing
    or malformed file — the converter then falls back to the parquet as-is.
    """
    refs: list[PsmRef] = []
    proteins: dict[str, str] = {}
    id_tags = {"PeptideIdentification", "UnassignedPeptideIdentification"}
    try:
        for _event, element in iterparse(consensusxml_path, events=("end",)):
            tag = _local(element.tag)
            if tag == "ProteinHit":
                pid = element.attrib.get("id")
                if pid:
                    proteins[pid] = element.attrib.get("accession", pid)
                element.clear()
            elif tag in id_tags:
                ref = _hit_ref(element, maplist, proteins)
                if ref is not None:
                    refs.append(ref)
                element.clear()
    except (OSError, ParseError, DefusedXmlException, KeyError, ValueError):
        return []
    return refs
