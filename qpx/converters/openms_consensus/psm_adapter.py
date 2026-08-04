"""consensusXML -> QPX psm records.

Each ``PeptideIdentification`` (assigned to a consensus feature or unassigned) is
one spectrum match. We emit one psm record per hit with PK
``[peptidoform, charge, run_file_name, scan]``. The run is resolved from the
identification's ``map_index`` / ``id_merge_index`` (→ the map's file) or, for a
single-run consensusXML, the sole run.
"""

from __future__ import annotations

import logging
import re

from qpx.converters.openms_consensus.feature_adapter import _run_stem, load_consensus_map, to_proforma

_log = logging.getLogger(__name__)

_SCAN_RE = re.compile(r"(?:scan|index|spectrum)=(\d+)", re.IGNORECASE)


def _scan_of(spectrum_ref: str) -> list[int]:
    """Parse the scan number(s) from a spectrum reference into a list<int>."""
    return [int(m) for m in _SCAN_RE.findall(str(spectrum_ref or ""))]


def _run_resolver(cm):
    """Build a callable mapping a PeptideIdentification to its run_file_name."""
    headers = cm.getColumnHeaders()
    map_run = {idx: _run_stem(headers[idx].filename) for idx in headers}
    distinct = sorted(set(map_run.values()))
    sole_run = distinct[0] if len(distinct) == 1 else None

    def resolve(pid) -> str | None:
        for key in ("map_index", "id_merge_index"):
            if pid.metaValueExists(key):
                run = map_run.get(int(pid.getMetaValue(key)))
                if run:
                    return run
        return sole_run

    return resolve


def _iter_peptide_ids(cm):
    """Yield every PeptideIdentification: unassigned + assigned to features."""
    yield from cm.getUnassignedPeptideIdentifications()
    for cf in cm:
        yield from cf.getPeptideIdentifications()


def consensus_psms_to_records(consensusxml_path: str | None = None, cm=None) -> list[dict]:
    """Return QPX psm record dicts extracted from a consensusXML.

    Pass either ``consensusxml_path`` (loaded here) or an already-loaded ``cm``.
    """
    cm = cm if cm is not None else load_consensus_map(consensusxml_path)
    resolve_run = _run_resolver(cm)

    records: list[dict] = []
    seen: set[tuple] = set()
    for pid in _iter_peptide_ids(cm):
        run = resolve_run(pid)
        if run is None:
            continue
        spectrum_ref = pid.getSpectrumReference() if hasattr(pid, "getSpectrumReference") else ""
        if not spectrum_ref and pid.metaValueExists("spectrum_reference"):
            spectrum_ref = pid.getMetaValue("spectrum_reference")
        scan = _scan_of(spectrum_ref)
        if not scan:
            # The PSM primary key is [peptidoform, charge, run_file_name, scan];
            # an identification whose spectrum_reference carries no scan token
            # cannot be keyed uniquely. Skip it rather than write a scan=[] record
            # that would collapse distinct spectra under the primary key.
            _log.debug("Skipping consensusXML PSM with no scan token in spectrum_reference: %r", spectrum_ref)
            continue
        obs_mz = float(pid.getMZ()) if pid.getMZ() else 0.0
        # When the identification score IS the q-value, the hit score is the
        # peptide q-value (OpenMS FDR output); otherwise it is a search score.
        score_type = str(pid.getScoreType() or "")
        score_is_qvalue = score_type.lower() in ("q-value", "qvalue", "fdr")
        for hit in pid.getHits():
            seq_obj = hit.getSequence()
            peptidoform = to_proforma(seq_obj)
            charge = int(hit.getCharge() or 0)
            calc_mz = float(seq_obj.getMZ(charge)) if charge else obs_mz
            key = (peptidoform, charge, run, tuple(scan))
            if key in seen:
                continue
            seen.add(key)
            is_decoy = hit.metaValueExists("target_decoy") and "decoy" in str(hit.getMetaValue("target_decoy")).lower()
            pep = None
            for mv in ("Posterior Error Probability_score", "PEP", "pep"):
                if hit.metaValueExists(mv):
                    pep = float(hit.getMetaValue(mv))
                    break
            score = float(hit.getScore()) if hit.getScore() is not None else None
            additional_scores = []
            if score is not None:
                # Route the identification score into additional_scores: the psm
                # schema has no dedicated q-value column.
                name = "q-value" if score_is_qvalue else (score_type or "search_score")
                additional_scores.append(
                    {"score_name": name, "score_value": score, "higher_better": bool(pid.isHigherScoreBetter())}
                )
            if hit.metaValueExists("consensus_support"):
                additional_scores.append(
                    {
                        "score_name": "consensus_support",
                        "score_value": float(hit.getMetaValue("consensus_support")),
                        "higher_better": True,
                    }
                )
            records.append(
                {
                    "sequence": seq_obj.toUnmodifiedString(),
                    "peptidoform": peptidoform,
                    "charge": charge,
                    "run_file_name": run,
                    "scan": scan,
                    "rt": float(pid.getRT()) if pid.getRT() else None,
                    "calculated_mz": calc_mz,
                    "observed_mz": obs_mz,
                    "posterior_error_probability": pep,
                    "additional_scores": additional_scores or None,
                    "is_decoy": is_decoy,
                }
            )
    return records
