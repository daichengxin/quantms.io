"""Orchestrate consensusXML (+ SDRF) -> a QPX dataset (feature/psm/pg [+ run/sample]).

Interim quantms path while OpenMS ``-out_qpx`` is pre-1.1. pg carries an interim
unnormalized unique-peptide-sum intensity (until ``-out_qpx`` provides the real
quant). run/sample come from the SDRF (reusing :class:`SdrfConverter`) when an
SDRF is provided; when it is, the consensusXML channels are also checked against
the SDRF ``comment[label]`` and mismatches are logged as warnings.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from qpx.converters.openms_consensus.feature_adapter import (
    check_channels_vs_sdrf,
    load_consensus_map,
)
from qpx.converters.openms_consensus.pg_adapter import accession_to_anchor, consensus_protein_groups_to_records
from qpx.core.data.identity import derive_id
from qpx.writers.feature import FeatureWriter
from qpx.writers.pg import PgWriter
from qpx.writers.psm import PsmWriter

_log = logging.getLogger(__name__)

_STRUCTURE_ALL = ("feature", "psm", "pg", "run", "sample")

# Producer-specific feature identity composite for the openms-consensus path
# (measured keys agreed in bigbio/qpx#229). Passed to the FeatureWriter so the
# derived feature_id hashes exactly these columns (all present in feature.yaml)
# instead of the schema default.
_FEATURE_IDENTITY_COMPOSITE = ("peptidoform", "charge", "run_file_name", "rt", "scan", "observed_mz")
# The psm view's schema identity_composite (psm.yaml) — used to derive the
# psm_id we cross-reference from feature.psm_ids / psm.feature_id.
_PSM_IDENTITY_COMPOSITE = ("peptidoform", "charge", "run_file_name", "scan")


def _link_feature_psm(feature_records: list[dict], psm_records: list[dict]) -> None:
    """Populate feature.psm_ids and psm.feature_id in place for one consensus feature.

    Ids are derived with the SAME composites the writers use, so these converter
    -populated cross-refs match the writer-derived primary keys byte-for-byte
    (``derive_id`` is deterministic). A consensus feature yields one feature
    record per run; each PSM links to the feature record of its own run. PSMs
    whose run has no feature record keep ``feature_id`` null (resolves #182).
    """
    feat_by_run: dict[str, tuple[int, list[int]]] = {}
    for rec in feature_records:
        feat_id = derive_id([rec.get(f) for f in _FEATURE_IDENTITY_COMPOSITE])
        feat_by_run[rec["run_file_name"]] = (feat_id, [])
    for prec in psm_records:
        entry = feat_by_run.get(prec.get("run_file_name"))
        if entry is None:
            continue
        feat_id, psm_ids = entry
        prec["feature_id"] = feat_id
        psm_id = derive_id([prec.get(f) for f in _PSM_IDENTITY_COMPOSITE])
        if psm_id not in psm_ids:
            psm_ids.append(psm_id)
    for rec in feature_records:
        _, psm_ids = feat_by_run[rec["run_file_name"]]
        if psm_ids:
            rec["psm_ids"] = psm_ids


# Above this consensusXML size, auto-select the streaming reader: the pyopenms
# in-memory load needs ~0.8x the file in RAM, so ~4 GB (≈3 GB RAM) is a safe point
# to switch to the low-memory path on a normal node.
_STREAM_THRESHOLD_BYTES = 4 * 1024**3


def _should_stream(path: str) -> bool:
    try:
        return Path(path).stat().st_size > _STREAM_THRESHOLD_BYTES
    except OSError:
        return False


def _convert_streaming(consensusxml_path, out, output_prefix, structures, sdrf_path, creator, pg_top) -> dict:
    """Single-pass, low-memory feature/psm/pg from a streamed consensusXML.

    One ordered ``iter_all()`` pass over the elements + unassigned IDs feeds the
    same per-element builders the in-memory adapters use (so output is identical),
    writing feature/psm in batches and accumulating the pg maps; pg records are
    built at the end. The whole map is never held, and the file is parsed once.
    """
    from collections import defaultdict
    from contextlib import ExitStack

    from qpx.converters.openms_consensus.feature_adapter import (
        _run_stem,
        feature_map_info,
        feature_records_for_cf,
    )
    from qpx.converters.openms_consensus.pg_adapter import (
        _ProteinMaps,
        accession_to_anchor,
        accumulate_cf_intensity,
        accumulate_cf_maps,
        accumulate_unassigned_maps,
        build_pg_records,
    )
    from qpx.converters.openms_consensus.psm_adapter import _run_resolver, psm_records_for_pid
    from qpx.converters.openms_consensus.streaming import StreamingConsensusMap

    cm = StreamingConsensusMap(consensusxml_path)
    if sdrf_path:
        for msg in check_channels_vs_sdrf(cm, sdrf_path):
            _log.warning("consensusXML/SDRF channel mismatch: %s", msg)

    map_info = feature_map_info(cm)
    headers = cm.getColumnHeaders()
    map_run = {i: _run_stem(headers[i].filename) for i in headers}
    want_feature, want_psm, want_pg = ("feature" in structures, "psm" in structures, "pg" in structures)
    anchor_map = accession_to_anchor(cm) if want_feature else None
    resolve_run = _run_resolver(cm) if want_psm else None
    maps = _ProteinMaps() if want_pg else None
    pep_intensity: dict = defaultdict(float) if want_pg else {}
    seen: set = set()
    written: dict[str, Path] = {}
    batch = 100_000

    with ExitStack() as stack:
        fw = pw = None
        if want_feature:
            fw = stack.enter_context(
                FeatureWriter(
                    str(out / f"{output_prefix}.feature.parquet"),
                    creator=creator,
                    identity_composite=_FEATURE_IDENTITY_COMPOSITE,
                )
            )
        if want_psm:
            pw = stack.enter_context(PsmWriter(str(out / f"{output_prefix}.psm.parquet"), creator=creator))
        feat_buf: list[dict] = []
        psm_buf: list[dict] = []
        for kind, obj in cm.iter_all():
            if kind == "element":
                cf_feats = feature_records_for_cf(obj, map_info, anchor_map) if fw is not None else []
                cf_psms: list[dict] = []
                if pw is not None:
                    for pid in obj.getPeptideIdentifications():
                        cf_psms.extend(psm_records_for_pid(pid, resolve_run, seen))
                # Cross-reference the feature<->psm ids only when both views are
                # emitted; each PSM links to the feature record of its own run.
                if fw is not None and pw is not None:
                    _link_feature_psm(cf_feats, cf_psms)
                feat_buf.extend(cf_feats)
                psm_buf.extend(cf_psms)
                if maps is not None:
                    accumulate_cf_maps(obj, map_run, maps)
                    accumulate_cf_intensity(obj, map_info, pep_intensity)
            else:  # unassigned peptide identification
                if pw is not None:
                    # Unassigned PSMs map to no feature -> feature_id stays null.
                    psm_buf.extend(psm_records_for_pid(obj, resolve_run, seen))
                if maps is not None:
                    accumulate_unassigned_maps(obj, map_run, maps)
            if fw is not None and len(feat_buf) >= batch:
                fw.write_batch(feat_buf)
                feat_buf = []
            if pw is not None and len(psm_buf) >= batch:
                pw.write_batch(psm_buf)
                psm_buf = []
        if fw is not None:
            if feat_buf:
                fw.write_batch(feat_buf)
            written["feature"] = out / f"{output_prefix}.feature.parquet"
        if pw is not None:
            if psm_buf:
                pw.write_batch(psm_buf)
            written["psm"] = out / f"{output_prefix}.psm.parquet"

    if want_pg:
        records = build_pg_records(cm, map_info, maps, pep_intensity, sdrf_path, pg_top)
        path = out / f"{output_prefix}.pg.parquet"
        with PgWriter(str(path), creator=creator) as w:
            if records:
                w.write_batch(records)
        written["pg"] = path
    return written


def _validate_structures(structures: tuple[str, ...], sdrf_path: Optional[str]) -> None:
    """Reject unknown structure names and run/sample requested without an SDRF."""
    unknown = [s for s in structures if s not in _STRUCTURE_ALL]
    if unknown:
        raise ValueError(f"Unknown structure(s) {unknown}; valid values are {list(_STRUCTURE_ALL)}")
    if ("run" in structures or "sample" in structures) and not sdrf_path:
        raise ValueError("The 'run'/'sample' structures require an SDRF (pass sdrf_path)")


class OpenMSConsensusConverter:  # pylint: disable=too-few-public-methods
    """consensusXML + SDRF -> QPX views.

    A single-entry orchestrator (``convert``) — the interim counterpart to the
    other converter classes, kept as a class for call-site symmetry with them.
    """

    def convert(
        self,
        consensusxml_path: str,
        output_folder: str,
        output_prefix: str = "openms",
        sdrf_path: Optional[str] = None,
        structures: tuple[str, ...] = _STRUCTURE_ALL,
        creator: str = "openms-consensus",
        pg_top: int = 0,
        streaming: Optional[bool] = None,
    ) -> dict[str, Path]:
        """Write the requested QPX views and return ``{structure: parquet path}``.

        ``structures`` selects which of feature/psm/pg/run/sample to emit. pg
        carries an interim unnormalized unique-peptide-sum intensity; ``pg_top``
        bounds the peptides used (0 = all; 3 mirrors ProteomicsLFQ/IsobaricWorkflow).

        ``streaming`` picks the consensusXML reader: ``None`` (default) auto-selects
        — the low-memory streaming reader for files above ``_STREAM_THRESHOLD_BYTES``
        (which pyopenms would otherwise load whole into ~0.8x-file RAM), else the
        faster in-memory pyopenms load. ``True``/``False`` forces the choice.
        """
        _validate_structures(structures, sdrf_path)

        out = Path(output_folder)
        out.mkdir(parents=True, exist_ok=True)
        written: dict[str, Path] = {}

        if {"feature", "psm", "pg"}.intersection(structures):
            use_stream = streaming if streaming is not None else _should_stream(consensusxml_path)
            if use_stream:
                # Low-memory path: a single ordered pass over the consensusXML builds
                # feature/psm/pg together and writes in batches, so the whole map is
                # never in memory and the file is parsed once (not once per view).
                written.update(_convert_streaming(consensusxml_path, out, output_prefix, structures, sdrf_path, creator, pg_top))
            else:
                # In-memory path: pyopenms loads the map once (fast for smaller files);
                # the adapters iterate it cheaply. Output is identical either way.
                written.update(
                    self._convert_in_memory(consensusxml_path, out, output_prefix, structures, sdrf_path, creator, pg_top)
                )

        if "run" in structures or "sample" in structures:
            from qpx.converters.sdrf import SdrfConverter

            with SdrfConverter() as sdrf_conv:
                sdrf_conv.convert(
                    sdrf_path=sdrf_path,
                    sample_output=str(out / f"{output_prefix}.sample.parquet"),
                    run_output=str(out / f"{output_prefix}.run.parquet"),
                )
            # Record only the structures the caller actually requested.
            if "run" in structures:
                written["run"] = out / f"{output_prefix}.run.parquet"
            if "sample" in structures:
                written["sample"] = out / f"{output_prefix}.sample.parquet"

        return written

    @staticmethod
    def _convert_in_memory(consensusxml_path, out, output_prefix, structures, sdrf_path, creator, pg_top) -> dict:
        """feature/psm/pg via the in-memory pyopenms map (loaded once, iterated cheaply)."""
        from qpx.converters.openms_consensus.feature_adapter import feature_map_info, feature_records_for_cf
        from qpx.converters.openms_consensus.psm_adapter import _run_resolver, psm_records_for_pid

        cm = load_consensus_map(consensusxml_path)
        if sdrf_path:
            for msg in check_channels_vs_sdrf(cm, sdrf_path):
                _log.warning("consensusXML/SDRF channel mismatch: %s", msg)
        written: dict[str, Path] = {}
        want_feature, want_psm = "feature" in structures, "psm" in structures
        if want_feature or want_psm:
            # Build feature and psm per consensus feature (mirroring the streaming
            # path's element loop) so their cross-refs can be linked identically:
            # assigned PSMs first (cf order), then the unassigned PSMs.
            map_info = feature_map_info(cm)
            # Share the protein-group leader map so feature.anchor_protein matches pg.
            anchor_map = accession_to_anchor(cm) if want_feature else None
            resolve_run = _run_resolver(cm) if want_psm else None
            seen: set = set()
            feat_recs: list[dict] = []
            psm_recs: list[dict] = []
            for cf in cm:
                cf_feats = feature_records_for_cf(cf, map_info, anchor_map) if want_feature else []
                cf_psms: list[dict] = []
                if want_psm:
                    for pid in cf.getPeptideIdentifications():
                        cf_psms.extend(psm_records_for_pid(pid, resolve_run, seen))
                if want_feature and want_psm:
                    _link_feature_psm(cf_feats, cf_psms)
                feat_recs.extend(cf_feats)
                psm_recs.extend(cf_psms)
            if want_psm:
                for pid in cm.getUnassignedPeptideIdentifications():
                    psm_recs.extend(psm_records_for_pid(pid, resolve_run, seen))
            if want_feature:
                path = out / f"{output_prefix}.feature.parquet"
                with FeatureWriter(str(path), creator=creator, identity_composite=_FEATURE_IDENTITY_COMPOSITE) as w:
                    if feat_recs:
                        w.write_batch(feat_recs)
                written["feature"] = path
            if want_psm:
                path = out / f"{output_prefix}.psm.parquet"
                with PsmWriter(str(path), creator=creator) as w:
                    if psm_recs:
                        w.write_batch(psm_recs)
                written["psm"] = path
        if "pg" in structures:
            recs = consensus_protein_groups_to_records(sdrf_path=sdrf_path, cm=cm, top=pg_top)
            path = out / f"{output_prefix}.pg.parquet"
            with PgWriter(str(path), creator=creator) as w:
                if recs:
                    w.write_batch(recs)
            written["pg"] = path
        return written
