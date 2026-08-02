"""consensusXML -> QPX feature records (core extraction).

Uses pyopenms to read the ConsensusMap. Each consensus feature is a peptide
linked across the map columns; a map column is one *run* (label-free) or one
*channel of a run* (isobaric, ``experiment_type == labeled_MS2``). We emit one
feature record per ``(peptidoform, charge, run_file_name, rt)`` carrying the
per-channel intensities as ``intensities: list<{label, intensity}>`` — matching
the QPX feature schema. Protein-level quantity is not produced here.
"""

from __future__ import annotations

import re
from pathlib import PurePosixPath
from typing import Optional

from qpx.converters.channel_labels import normalize_label

# OpenMS isobaric channel labels look like ``tmt6plex_126`` / ``itraq4plex_114``.
_CHANNEL_RE = re.compile(r"(tmt|itraq)\d*plex?_(\d+[NC]?)", re.IGNORECASE)


def _canonical_channel(label: Optional[str]) -> str:
    """Canonicalize an OpenMS map label to the QPX channel label.

    ``tmt6plex_126`` -> ``TMT126``; ``itraq4plex_114`` -> ``ITRAQ114``. Falls back
    to :func:`normalize_label`, then the raw string, so unknown labels are kept
    rather than dropped.
    """
    if not label:
        return "LFQ"
    m = _CHANNEL_RE.search(str(label))
    if m:
        family = "TMT" if m.group(1).lower() == "tmt" else "ITRAQ"
        return normalize_label(f"{family}{m.group(2).upper()}")
    return normalize_label(str(label))


def _run_stem(filename: str) -> str:
    """Run file name = basename with the final extension stripped."""
    name = PurePosixPath(str(filename).replace("\\", "/")).name
    return re.sub(r"(?i)\.(mzml|mzml\.gz|raw|d|wiff|mgf)$", "", name)


def to_proforma(aa_sequence) -> str:
    """Convert an OpenMS ``AASequence`` to a ProForma peptidoform string.

    OpenMS ``toUniModString()`` uses ``(UniMod:N)`` with a leading/trailing ``.``
    for terminal mods (e.g. ``.(UniMod:737)THSQEEM(UniMod:35)QHMQR``); ProForma
    uses ``[UNIMOD:N]`` and ``[..]-`` / ``-[..]`` for N-/C-terminal mods
    (``[UNIMOD:737]-THSQEEM[UNIMOD:35]QHMQR``).
    """
    s = aa_sequence.toUniModString().replace("UniMod:", "UNIMOD:")
    s = re.sub(r"^\.\(([^)]+)\)", r"[\1]-", s)  # N-terminal
    s = re.sub(r"\.\(([^)]+)\)$", r"-[\1]", s)  # C-terminal
    s = re.sub(r"\(([^)]+)\)", r"[\1]", s)  # internal residue mods
    return s


def load_consensus_map(consensusxml_path: str):
    """Parse a consensusXML into a ``ConsensusMap`` — the single parse point.

    The feature/psm/pg adapters each only *read* the map, so one loaded instance
    is shared across all three (see ``OpenMSConsensusConverter.convert``) instead
    of re-parsing the (possibly large) XML three times.
    """
    import pyopenms as oms

    cm = oms.ConsensusMap()
    oms.ConsensusXMLFile().load(str(consensusxml_path), cm)
    return cm


def consensus_features_to_records(consensusxml_path: str | None = None, cm=None) -> list[dict]:
    """Return QPX feature record dicts extracted from a consensusXML.

    Pass either ``consensusxml_path`` (loaded here) or an already-loaded ``cm``.
    """
    cm = cm if cm is not None else load_consensus_map(consensusxml_path)

    headers = cm.getColumnHeaders()
    is_labeled = cm.getExperimentType() != "label_free"
    map_info: dict[int, tuple[str, str]] = {}
    for idx in headers:
        header = headers[idx]
        run = _run_stem(header.filename)
        label = _canonical_channel(header.label) if is_labeled else "LFQ"
        map_info[idx] = (run, label)

    records: list[dict] = []
    for cf in cm:
        pids = cf.getPeptideIdentifications()
        if not pids or not pids[0].getHits():
            continue
        spectrum_ref = pids[0].getSpectrumReference() if hasattr(pids[0], "getSpectrumReference") else ""
        if not spectrum_ref and pids[0].metaValueExists("spectrum_reference"):
            spectrum_ref = pids[0].getMetaValue("spectrum_reference")
        scan = [int(m) for m in re.findall(r"(?:scan|index|spectrum)=(\d+)", str(spectrum_ref or ""), re.IGNORECASE)]
        hit = pids[0].getHits()[0]
        seq_obj = hit.getSequence()
        peptidoform = to_proforma(seq_obj)
        sequence = seq_obj.toUnmodifiedString()
        charge = int(cf.getCharge() or hit.getCharge() or 0)
        is_decoy = False
        if hit.metaValueExists("target_decoy"):
            is_decoy = "decoy" in str(hit.getMetaValue("target_decoy")).lower()
        observed_mz = float(cf.getMZ()) if cf.getMZ() else 0.0
        calculated_mz = float(seq_obj.getMZ(charge)) if charge else observed_mz
        evidences = hit.getPeptideEvidences()
        anchor_protein = evidences[0].getProteinAccession() if evidences else None
        if isinstance(anchor_protein, bytes):
            anchor_protein = anchor_protein.decode()

        # Group the sub-feature intensities by run: for isobaric, one run carries
        # several channel maps (same rt); for label-free, one map == one run. One
        # entry per label — if a channel is linked twice, keep the max rather than
        # emit the channel twice (which would double-count downstream).
        by_run: dict[str, dict] = {}
        for sub in cf.getFeatureList():
            intensity = float(sub.getIntensity())
            if intensity <= 0:
                continue
            run, label = map_info.get(sub.getMapIndex(), (None, None))
            if run is None:
                continue
            entry = by_run.setdefault(run, {"rt": float(sub.getRT()), "labels": {}})
            entry["labels"][label] = max(entry["labels"].get(label, 0.0), intensity)

        for run, entry in by_run.items():
            intensities = [{"label": label, "intensity": inten} for label, inten in entry["labels"].items()]
            records.append(
                {
                    "sequence": sequence,
                    "peptidoform": peptidoform,
                    "charge": charge,
                    "run_file_name": run,
                    "rt": entry["rt"],
                    "scan": scan,
                    "intensities": intensities,
                    "is_decoy": is_decoy,
                    "calculated_mz": calculated_mz,
                    "observed_mz": observed_mz,
                    "anchor_protein": anchor_protein,
                }
            )
    return records
