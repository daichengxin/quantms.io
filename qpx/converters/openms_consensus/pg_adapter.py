"""consensusXML -> QPX pg (protein group) records — identification only.

Protein groups come from the OpenMS protein-inference graph: indistinguishable
protein groups when present, else each protein hit as a singleton group. Peptide
counts come from the peptide→protein evidence links.

**No protein intensity is emitted** (interim): the consensusXML has no
protein-level abundance. Rows carry ``label``/``intensity`` null (identification-
only) until OpenMS ``-out_qpx`` provides the authoritative protein quant. One row
per ``(protein group, grouped_runs)`` — the quantification unit(s), from the SDRF
when given, otherwise the whole run set as a single unit.
"""

from __future__ import annotations

import re
from collections import defaultdict
from typing import Optional

from qpx.converters.channel_labels import fraction_groups_from_sdrf
from qpx.converters.openms_consensus.feature_adapter import (
    _canonical_channel,
    _run_stem,
    load_consensus_map,
    to_proforma,
)

_GENE_RE = re.compile(r"GN=([^\s]+)")


def _protein_maps(cm) -> tuple[dict[str, set[str]], dict[str, set[str]], dict[str, set[tuple]]]:
    """Return (accession -> peptide sequences, -> runs identified in, -> features).

    The per-protein run set lets pg rows be emitted only for the quantification
    units where the protein was actually identified; the feature set (distinct
    peptidoform+charge) feeds ``feature_counts``.
    """
    headers = cm.getColumnHeaders()
    map_run = {i: _run_stem(headers[i].filename) for i in headers}
    acc_to_pep: dict[str, set[str]] = defaultdict(set)
    acc_to_runs: dict[str, set[str]] = defaultdict(set)
    acc_to_feat: dict[str, set[tuple]] = defaultdict(set)

    def _accession(ev) -> str | None:
        acc = ev.getProteinAccession()
        if not acc:
            return None
        return acc.decode() if isinstance(acc, bytes) else str(acc)

    def _collect(pid, runs):
        for hit in pid.getHits():
            seq_obj = hit.getSequence()
            seq = seq_obj.toUnmodifiedString()
            feat = (to_proforma(seq_obj), int(hit.getCharge() or 0))
            for ev in hit.getPeptideEvidences():
                acc = _accession(ev)
                if acc:
                    acc_to_pep[acc].add(seq)
                    acc_to_feat[acc].add(feat)
                    acc_to_runs[acc] |= runs

    # Assigned IDs: the runs are the consensus feature's member maps.
    for cf in cm:
        cf_runs = {map_run.get(sub.getMapIndex()) for sub in cf.getFeatureList() if sub.getIntensity() > 0}
        cf_runs.discard(None)
        for pid in cf.getPeptideIdentifications():
            _collect(pid, cf_runs)
    # Unassigned IDs: the run comes from map_index / id_merge_index.
    for pid in cm.getUnassignedPeptideIdentifications():
        run = None
        for key in ("map_index", "id_merge_index"):
            if pid.metaValueExists(key):
                run = map_run.get(int(pid.getMetaValue(key)))
                break
        _collect(pid, {run} if run else set())
    return acc_to_pep, acc_to_runs, acc_to_feat


def _is_decoy_accession(acc: str) -> bool:
    return str(acc).upper().startswith(("DECOY", "REV_", "RANDOM_"))


def _is_contaminant(acc: str) -> bool:
    return "CONTAM" in str(acc).upper()


def _acc_str(acc) -> str:
    return acc.decode() if isinstance(acc, bytes) else str(acc)


def _protein_hit_meta(prot) -> tuple[dict[str, bool], dict[str, float], dict[str, str]]:
    """Per-accession decoy flag, q-value (when the protein score IS a q-value), gene."""
    score_is_qvalue = str(prot.getScoreType() or "").lower() in ("q-value", "qvalue", "fdr")
    acc_decoy: dict[str, bool] = {}
    acc_qvalue: dict[str, float] = {}
    acc_gene: dict[str, str] = {}
    for hit in prot.getHits():
        acc = _acc_str(hit.getAccession())
        if hit.metaValueExists("target_decoy"):
            acc_decoy[acc] = "decoy" in str(hit.getMetaValue("target_decoy")).lower()
        if score_is_qvalue and hit.getScore() is not None:
            acc_qvalue[acc] = float(hit.getScore())
        gene = _GENE_RE.search(str(hit.getDescription() or "")) if hasattr(hit, "getDescription") else None
        if gene:
            acc_gene[acc] = gene.group(1)
    return acc_decoy, acc_qvalue, acc_gene


def _build_groups(prot) -> list[list[str]]:
    """Protein groups: indistinguishable groups if present, else one per hit."""
    indistinguishable = prot.getIndistinguishableProteins()
    if indistinguishable:
        return [[_acc_str(a) for a in grp.accessions] for grp in indistinguishable if grp.accessions]
    return [[_acc_str(hit.getAccession())] for hit in prot.getHits()]


def consensus_protein_groups_to_records(
    consensusxml_path: str | None = None,
    sdrf_path: Optional[str] = None,
    cm=None,
) -> list[dict]:
    """Return QPX pg record dicts: one per (protein group, unit, label), null intensity.

    The label (channel/sample) dimension is populated so each quantification slot
    exists with a placeholder null ``intensity`` — filled later by the OpenMS-team
    ``-out_qpx``. No protein quant is invented here. Pass either
    ``consensusxml_path`` (loaded here) or an already-loaded ``cm``.
    """
    cm = cm if cm is not None else load_consensus_map(consensusxml_path)

    headers = cm.getColumnHeaders()
    all_runs = sorted({_run_stem(headers[i].filename) for i in headers})
    # grouped_runs unit per run, from the SDRF (fractions grouped); else one unit.
    run_to_grouped = fraction_groups_from_sdrf(sdrf_path)
    if run_to_grouped:
        units = {tuple(v) for v in run_to_grouped.values()}
    else:
        units = {tuple(all_runs)}

    # The label (channel/sample) dimension: one pg row per protein x unit x label
    # so the quantification slot exists with a null intensity — a placeholder the
    # OpenMS-team -out_qpx will later fill. Labels come from the consensusXML map
    # columns (isobaric channels), or "LFQ" for label-free.
    if cm.getExperimentType() != "label_free":
        labels = sorted({_canonical_channel(headers[i].label) for i in headers})
    else:
        labels = ["LFQ"]

    acc_to_pep, acc_to_runs, acc_to_feat = _protein_maps(cm)

    # Protein groups + per-accession decoy/q-value/gene from the inference graph.
    prot_ids = cm.getProteinIdentifications()
    if prot_ids:
        acc_decoy, acc_qvalue, acc_gene = _protein_hit_meta(prot_ids[0])
        groups = _build_groups(prot_ids[0])
    else:
        acc_decoy, acc_qvalue, acc_gene, groups = {}, {}, {}, []

    records: list[dict] = []
    for accs in groups:
        anchor = accs[0]
        peptide_seqs: set[str] = set()
        feats: set[tuple] = set()
        for acc in accs:
            peptide_seqs |= acc_to_pep.get(acc, set())
            feats |= acc_to_feat.get(acc, set())
        n_pep = len(peptide_seqs)
        n_feat = len(feats)
        # Prefer the target_decoy meta; fall back to the accession prefix.
        is_decoy = all(acc_decoy.get(a, _is_decoy_accession(a)) for a in accs)
        qvals = [acc_qvalue[a] for a in accs if a in acc_qvalue]
        global_qvalue = min(qvals) if qvals else None
        genes = [acc_gene[a] for a in accs if a in acc_gene] or None
        # Only the quantification units where this group was actually identified
        # (its peptides appear in a run of that unit) — not every unit.
        group_runs: set[str] = set()
        for acc in accs:
            group_runs |= acc_to_runs.get(acc, set())
        group_units = [unit for unit in units if group_runs.intersection(unit)] or list(units)
        for unit in group_units:
            for label in labels:
                records.append(
                    {
                        "pg_accessions": list(accs),
                        "anchor_protein": anchor,
                        "grouped_runs": list(unit),
                        "label": label,  # channel/sample slot; intensity filled later
                        "intensity": None,  # interim: no protein quant yet (OpenMS -out_qpx)
                        "global_qvalue": global_qvalue,
                        "is_decoy": is_decoy,
                        "contaminant": any(_is_contaminant(a) for a in accs),
                        "gg_accessions": genes,
                        "gg_names": genes,
                        "peptide_counts": {"unique_sequences": n_pep, "total_sequences": n_pep},
                        "feature_counts": {"unique_features": n_feat, "total_features": n_feat},
                        "peptides": [{"protein_name": a, "peptide_count": n_pep} for a in accs],
                    }
                )
    return records
