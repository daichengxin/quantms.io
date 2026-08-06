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

import pyarrow as pa

from qpx.converters.openms_consensus.feature_adapter import (
    check_channels_vs_sdrf,
    load_consensus_map,
)
from qpx.converters.openms_consensus.pg_adapter import accession_to_anchor, consensus_protein_groups_to_records
from qpx._version import __version__
from qpx.converters.orchestrator import BaseOrchestrator
from qpx.core.constants import FEATURE, ONTOLOGY, PG, PSM, RUN, SAMPLE
from qpx.core.data import FeatureSchema, PsmSchema
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
_FEATURE_IDENTITY_COMPOSITE = (
    "peptidoform",
    "charge",
    "run_file_name",
    "rt",
    "scan",
    "observed_mz",
    "consensus_rt",
)
# The psm view's schema identity_composite (psm.yaml) — used to derive the
# psm_id we cross-reference from feature.psm_ids / psm.feature_id.
_PSM_IDENTITY_COMPOSITE = ("peptidoform", "charge", "run_file_name", "scan")


def _derive_persisted_record_id(record: dict, composite: tuple[str, ...], schema: pa.Schema) -> int:
    """Derive an id from the values exactly as they will be stored by Arrow."""
    values = [pa.scalar(record.get(name), type=schema.field(name).type).as_py() for name in composite]
    return derive_id(values)


def _link_feature_psm(feature_records: list[dict], psm_records: list[dict]) -> None:
    """Populate feature.psm_ids and psm.feature_id in place for one consensus feature.

    Ids are derived with the SAME composites the writers use, so these converter
    -populated cross-refs match the writer-derived primary keys byte-for-byte
    (``derive_id`` is deterministic). A consensus feature yields one feature
    record per run; each PSM links to the feature record of its own run. PSMs
    whose run has no feature record keep ``feature_id`` null (resolves #182).
    """
    feature_schema = FeatureSchema.get_arrow_schema()
    psm_schema = PsmSchema.get_arrow_schema()
    feat_by_run: dict[str, tuple[int, list[int]]] = {}
    for rec in feature_records:
        feat_id = _derive_persisted_record_id(rec, _FEATURE_IDENTITY_COMPOSITE, feature_schema)
        feat_by_run[rec["run_file_name"]] = (feat_id, [])
    for prec in psm_records:
        entry = feat_by_run.get(prec.get("run_file_name"))
        if entry is None:
            continue
        feat_id, psm_ids = entry
        prec["feature_id"] = feat_id
        psm_id = _derive_persisted_record_id(prec, _PSM_IDENTITY_COMPOSITE, psm_schema)
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


def _convert_streaming(consensusxml_path, out, output_prefix, structures, sdrf_path, creator, pg_top, compression="zstd") -> dict:
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
    from qpx.converters.openms_consensus.psm_adapter import _cf_element_runs, _run_resolver, psm_records_for_pid
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
    resolve_run = _run_resolver(cm) if want_psm or want_pg else None
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
                    compression=compression,
                )
            )
        if want_psm:
            pw = stack.enter_context(
                PsmWriter(
                    str(out / f"{output_prefix}.psm.parquet"),
                    creator=creator,
                    identity_composite=_PSM_IDENTITY_COMPOSITE,
                    compression=compression,
                )
            )
        feat_buf: list[dict] = []
        psm_buf: list[dict] = []
        for kind, obj in cm.iter_all():
            if kind == "element":
                cf_feats = feature_records_for_cf(obj, map_info, anchor_map) if fw is not None else []
                cf_psms: list[dict] = []
                if pw is not None:
                    # Multi-run isobaric PIDs carry a local id_merge_index; the
                    # feature's element runs disambiguate which run they belong to.
                    cf_runs = _cf_element_runs(obj, map_info)
                    for pid in obj.getPeptideIdentifications():
                        cf_psms.extend(psm_records_for_pid(pid, resolve_run, seen, cf_runs=cf_runs))
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
                    accumulate_unassigned_maps(obj, resolve_run, maps)
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
        with PgWriter(str(path), creator=creator, compression=compression) as w:
            if records:
                w.write_batch(records)
        written["pg"] = path
    return written


def _validate_structures(structures: tuple[str, ...], sdrf_path: Optional[str]) -> None:
    """Reject unknown structure names and run/sample requested without an SDRF."""
    unknown = [s for s in structures if s not in _STRUCTURE_ALL]
    if unknown:
        raise ValueError(f"Unknown structure(s) {unknown}; valid values are {list(_STRUCTURE_ALL)}")
    metadata = set(structures).intersection({"run", "sample"})
    if metadata and not sdrf_path:
        raise ValueError(f"An SDRF is required to write {sorted(metadata)}")


def _collect_score_names(table_path: Path) -> set[str]:
    """Read ``additional_scores`` score names from a written parquet file."""
    table = pa.parquet.read_table(str(table_path), columns=["additional_scores"])
    names: set[str] = set()
    for row in table.column("additional_scores").to_pylist():
        if row:
            for score in row:
                if score and score.get("score_name"):
                    names.add(score["score_name"])
    return names


def _consensus_is_isobaric(consensusxml_path: str) -> bool:
    """Return whether a consensusXML's map labels indicate an isobaric (TMT/iTRAQ) run."""
    try:
        cm = load_consensus_map(consensusxml_path)
        headers = cm.getColumnHeaders()
        from qpx.converters.openms_consensus.feature_adapter import consensus_channels

        return bool(consensus_channels(cm))
    except Exception:  # noqa: BLE001 - best-effort provenance hint
        return False


def _write_sdrf_metadata(
    output_folder: Path,
    output_prefix: str,
    sdrf_path: str,
    requested: set[str],
    compression: str = "zstd",
) -> tuple[dict[str, Path], list[dict]]:
    """Write the requested SDRF-backed run/sample structures.

    Returns ``(paths, run_ontology_entries)`` — the latter collected from the
    same ``SdrfConverter`` instance (its parse state populates the run ontology),
    so callers need not re-parse the SDRF.
    """
    metadata = requested.intersection({"run", "sample"})
    if not metadata:
        return {}, []

    from qpx.converters.sdrf import SdrfConverter

    paths = {
        "sample": output_folder / f"{output_prefix}.sample.parquet",
        "run": output_folder / f"{output_prefix}.run.parquet",
    }
    with SdrfConverter(compression=compression) as sdrf_converter:
        sdrf_converter.convert(
            sdrf_path=sdrf_path,
            sample_output=str(paths["sample"]) if "sample" in metadata else None,
            run_output=str(paths["run"]) if "run" in metadata else None,
        )
        run_ontology = list(sdrf_converter.run_ontology_entries())
    return {name: paths[name] for name in metadata}, run_ontology


class OpenMSConsensusConverter(BaseOrchestrator):  # pylint: disable=too-few-public-methods
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
        project_accession: Optional[str] = None,
        compression: str = "zstd",
    ) -> dict[str, Path]:
        """Write the requested QPX views and return ``{structure: parquet path}``.

        ``structures`` selects which of feature/psm/pg/run/sample to emit. pg
        carries an interim unnormalized unique-peptide-sum intensity; ``pg_top``
        bounds the peptides used (0 = all; 3 mirrors ProteomicsLFQ/IsobaricWorkflow).

        ``streaming`` picks the consensusXML reader: ``None`` (default) auto-selects
        — the low-memory streaming reader for files above ``_STREAM_THRESHOLD_BYTES``
        (which pyopenms would otherwise load whole into ~0.8x-file RAM), else the
        faster in-memory pyopenms load. ``True``/``False`` forces the choice.

        ``project_accession`` (e.g. ``PXD001819``) is stamped into dataset.parquet.

        ``compression`` is the Parquet codec for every written view (default zstd).
        """
        self._compression = compression
        _validate_structures(structures, sdrf_path)
        requested = set(structures)

        out = Path(output_folder)
        out.mkdir(parents=True, exist_ok=True)
        written: dict[str, Path] = {}

        if {"feature", "psm", "pg"}.intersection(structures):
            use_stream = streaming if streaming is not None else _should_stream(consensusxml_path)
            if use_stream:
                # Low-memory path: a single ordered pass over the consensusXML builds
                # feature/psm/pg together and writes in batches, so the whole map is
                # never in memory and the file is parsed once (not once per view).
                written.update(_convert_streaming(consensusxml_path, out, output_prefix, structures, sdrf_path, creator, pg_top, compression))
            else:
                # In-memory path: pyopenms loads the map once (fast for smaller files);
                # the adapters iterate it cheaply. Output is identical either way.
                written.update(
                    self._convert_in_memory(consensusxml_path, out, output_prefix, structures, sdrf_path, creator, pg_top, compression)
                )

        sdrf_paths, run_ontology = _write_sdrf_metadata(out, output_prefix, sdrf_path, requested, compression)
        written.update(sdrf_paths)

        # Metadata tables (ontology / provenance / dataset) — written when the
        # core + run/sample structures they describe are present. Best-effort:
        # ontology entries only materialise when a PSI-MS term resolves.
        if written:
            self._write_metadata(out, output_prefix, written, consensusxml_path, run_ontology, requested, project_accession)

        return written

    def _write_metadata(
        self,
        out: Path,
        output_prefix: str,
        written: dict[str, Path],
        consensusxml_path: str,
        run_ontology: list[dict],
        requested: set[str],
        project_accession: Optional[str] = None,
    ) -> None:
        """Write ontology/provenance/dataset metadata tables (best-effort).

        Mirrors the ``openms`` enrichment path: ontology collects the run-level
        terms from the SDRF plus the score names discovered in the written core
        tables; provenance records the consensusXML -> QPX conversion steps; the
        dataset table carries the project-level metadata.
        """
        ontology_entries = list(run_ontology)
        for view in (PSM, FEATURE, PG):
            path = written.get(view)
            if not path:
                continue
            try:
                names = _collect_score_names(path)
            except (KeyError, pa.ArrowInvalid):
                continue
            if names:
                from qpx.core.scores import score_ontology_entries

                ontology_entries.extend(score_ontology_entries(names, view=view))  # noqa: PERF401
        self._write_ontology(out, output_prefix, ontology_entries)

        structures = sorted(requested & {PSM, FEATURE, PG, RUN, SAMPLE})
        provenance_records = self._build_provenance(structures, consensusxml_path)
        self._write_provenance(out, output_prefix, provenance_records)
        if provenance_records:
            self._write_dataset(
                out,
                output_prefix,
                project_accession,
                software_name=provenance_records[0]["tool_name"],
                software_version=None,
                provenance_records=provenance_records,
            )

    @staticmethod
    def _build_provenance(structures: list[str], consensusxml_path: str) -> list[dict]:
        """Provenance records: consensusXML parsing + QPX conversion."""
        is_isobaric = _consensus_is_isobaric(consensusxml_path)
        step_name = "isobaric_quantification" if is_isobaric else "label_free_quantification"
        tool_name = "OpenMS/IsobaricWorkflow" if is_isobaric else "OpenMS/ProteomicsLFQ"
        return [
            {
                "step_order": 1,
                "step_category": "quantification",
                "step_name": step_name,
                "tool_name": tool_name,
                "tool_version": None,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                "output_views": [v for v in (FEATURE, PSM, PG) if v in structures],
            },
            {
                "step_order": 2,
                "step_category": "format_conversion",
                "step_name": "openms_consensus_qpx_conversion",
                "tool_name": "qpx",
                "tool_version": __version__,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                "output_views": [SAMPLE, RUN, ONTOLOGY],
            },
        ]

    @staticmethod
    def _convert_in_memory(consensusxml_path, out, output_prefix, structures, sdrf_path, creator, pg_top, compression="zstd") -> dict:
        """feature/psm/pg via the in-memory pyopenms map (loaded once, iterated cheaply)."""
        from qpx.converters.openms_consensus.feature_adapter import feature_map_info, feature_records_for_cf
        from qpx.converters.openms_consensus.psm_adapter import _cf_element_runs, _run_resolver, psm_records_for_pid

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
                    # Multi-run isobaric PIDs carry a local id_merge_index; the
                    # feature's element runs disambiguate which run they belong to.
                    cf_runs = _cf_element_runs(cf, map_info)
                    for pid in cf.getPeptideIdentifications():
                        cf_psms.extend(psm_records_for_pid(pid, resolve_run, seen, cf_runs=cf_runs))
                if want_feature and want_psm:
                    _link_feature_psm(cf_feats, cf_psms)
                feat_recs.extend(cf_feats)
                psm_recs.extend(cf_psms)
            if want_psm:
                for pid in cm.getUnassignedPeptideIdentifications():
                    psm_recs.extend(psm_records_for_pid(pid, resolve_run, seen))
            if want_feature:
                path = out / f"{output_prefix}.feature.parquet"
                with FeatureWriter(str(path), creator=creator, identity_composite=_FEATURE_IDENTITY_COMPOSITE, compression=compression) as w:
                    if feat_recs:
                        w.write_batch(feat_recs)
                written["feature"] = path
            if want_psm:
                path = out / f"{output_prefix}.psm.parquet"
                with PsmWriter(str(path), creator=creator, identity_composite=_PSM_IDENTITY_COMPOSITE, compression=compression) as w:
                    if psm_recs:
                        w.write_batch(psm_recs)
                written["psm"] = path
        if "pg" in structures:
            recs = consensus_protein_groups_to_records(sdrf_path=sdrf_path, cm=cm, top=pg_top)
            path = out / f"{output_prefix}.pg.parquet"
            with PgWriter(str(path), creator=creator, compression=compression) as w:
                if recs:
                    w.write_batch(recs)
            written["pg"] = path
        return written
