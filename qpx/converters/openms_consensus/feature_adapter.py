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


def _map_label(raw_label) -> str:
    """Resolve a consensusXML map label to a QPX channel, or ``LFQ``.

    Detection is by the label itself, NOT ``experiment_type``: quantms/OpenMS
    write ``experiment_type="label-free"`` even for isobaric runs (the maps still
    carry ``tmt6plex_126`` / ``itraq4plex_114`` labels; label-free maps carry the
    literal ``"label-free"``). Only a real isobaric channel pattern becomes a
    channel; everything else (``"label-free"``, empty, unknown) is ``LFQ``.
    """
    return _canonical_channel(raw_label) if _CHANNEL_RE.search(str(raw_label or "")) else "LFQ"


def consensus_channels(cm) -> set[str]:
    """Canonical isobaric channels present in the consensusXML maps (excludes LFQ)."""
    headers = cm.getColumnHeaders()
    return {lab for i in headers if (lab := _map_label(headers[i].label)) != "LFQ"}


def check_channels_vs_sdrf(cm, sdrf_path) -> list[str]:
    """Return warnings where the consensusXML channels disagree with the SDRF.

    The consensusXML map labels are the channel identity the converter uses; the
    SDRF ``comment[label]`` set is the experiment's declared ground truth. A
    channel present in one but not the other is a metadata inconsistency worth
    surfacing (e.g. a mis-declared plex, or the wrong SDRF). Returns human-readable
    messages — empty when the two agree, or when the SDRF has no label column.
    """
    from qpx.converters.channel_labels import normalize_label, read_sdrf_labels

    sdrf_raw = read_sdrf_labels(sdrf_path)
    if not sdrf_raw:
        return []
    sdrf_channels = {normalize_label(v) for v in sdrf_raw} - {"LFQ"}
    cm_channels = consensus_channels(cm)
    messages: list[str] = []
    only_cm = sorted(cm_channels - sdrf_channels)
    only_sdrf = sorted(sdrf_channels - cm_channels)
    if only_cm:
        messages.append(f"channel(s) in the consensusXML but not the SDRF comment[label]: {only_cm}")
    if only_sdrf:
        messages.append(f"channel(s) in the SDRF comment[label] but not the consensusXML: {only_sdrf}")
    return messages


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


def _pid_scans(pid) -> list[int]:
    """Scan numbers parsed from one identification's spectrum reference."""
    ref = pid.getSpectrumReference() if hasattr(pid, "getSpectrumReference") else ""
    if not ref and pid.metaValueExists("spectrum_reference"):
        ref = pid.getMetaValue("spectrum_reference")
    return [int(m) for m in re.findall(r"(?:scan|index|spectrum)=(\d+)", str(ref or ""), re.IGNORECASE)]


def _scan_by_run(pids, map_info: dict[int, tuple[str, str]]) -> dict[str, list[int]]:
    """Attribute each identification's scan(s) to its own run.

    A consensus feature links spectra from several runs, so scans are resolved
    per ID via its ``map_index`` (falling back to the sole run only when every
    map is the same physical run, e.g. isobaric channels) rather than copying
    one ID's scan onto every run's record.
    """
    runs_in_map = {info[0] for info in map_info.values()}
    sole_run = next(iter(runs_in_map)) if len(runs_in_map) == 1 else None
    scan_by_run: dict[str, list[int]] = {}
    for pid in pids:
        scans = _pid_scans(pid)
        if not scans:
            continue
        pid_run = None
        if pid.metaValueExists("map_index"):
            pid_run = map_info.get(int(pid.getMetaValue("map_index")), (None, None))[0]
        pid_run = pid_run or sole_run
        if pid_run is not None:
            scan_by_run.setdefault(pid_run, scans)
    return scan_by_run


def _group_subfeatures_by_run(cf, map_info: dict[int, tuple[str, str]]) -> dict[str, dict]:
    """Group a consensus feature's sub-feature intensities by run.

    For isobaric, one run carries several channel maps (same rt); for label-free,
    one map == one run. One entry per label — if a channel is linked twice, keep
    the max rather than emit it twice (which would double-count downstream).
    """
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
    return by_run


def consensus_features_to_records(consensusxml_path: str | None = None, cm=None, anchor_map=None) -> list[dict]:
    """Return QPX feature record dicts extracted from a consensusXML.

    Pass either ``consensusxml_path`` (loaded here) or an already-loaded ``cm``.
    ``anchor_map`` (accession -> protein-group leader, from
    ``pg_adapter.accession_to_anchor``) makes each feature's ``anchor_protein``
    the same group leader the pg view uses, so the feature->pg join is reliable;
    without it the anchor falls back to the peptide's first protein evidence.
    """
    cm = cm if cm is not None else load_consensus_map(consensusxml_path)
    map_info = feature_map_info(cm)
    records: list[dict] = []
    for cf in cm:
        records.extend(feature_records_for_cf(cf, map_info, anchor_map))
    return records


def feature_map_info(cm) -> dict[int, tuple[str, str]]:
    """Map index -> (run_file_name, label) from the consensusXML column headers.

    Label comes from the map label itself (isobaric channel -> canonical, else
    LFQ) — experiment_type is unreliable (quantms writes "label-free" for TMT).
    """
    headers = cm.getColumnHeaders()
    return {idx: (_run_stem(headers[idx].filename), _map_label(headers[idx].label)) for idx in headers}


def feature_records_for_cf(cf, map_info: dict[int, tuple[str, str]], anchor_map=None) -> list[dict]:
    """Feature records for one consensus feature (one per run, channels as intensities)."""
    pids = cf.getPeptideIdentifications()
    if not pids or not pids[0].getHits():
        return []
    scan_by_run = _scan_by_run(pids, map_info)
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
    # Resolve to the protein-group leader so feature.anchor_protein matches the
    # pg view (peptide evidence order alone does not identify the leader).
    if anchor_map and anchor_protein is not None:
        anchor_protein = anchor_map.get(anchor_protein, anchor_protein)

    records: list[dict] = []
    for run, entry in _group_subfeatures_by_run(cf, map_info).items():
        intensities = [{"label": label, "intensity": inten} for label, inten in entry["labels"].items()]
        records.append(
            {
                "sequence": sequence,
                "peptidoform": peptidoform,
                "charge": charge,
                "run_file_name": run,
                "scan": scan_by_run.get(run, []),
                "rt": entry["rt"],
                "intensities": intensities,
                "is_decoy": is_decoy,
                "calculated_mz": calculated_mz,
                "observed_mz": observed_mz,
                "anchor_protein": anchor_protein,
            }
        )
    return records
