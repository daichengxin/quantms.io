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
    consensus_features_to_records,
    load_consensus_map,
)
from qpx.converters.openms_consensus.pg_adapter import accession_to_anchor, consensus_protein_groups_to_records
from qpx.converters.openms_consensus.psm_adapter import consensus_psms_to_records
from qpx.writers.feature import FeatureWriter
from qpx.writers.pg import PgWriter
from qpx.writers.psm import PsmWriter

_log = logging.getLogger(__name__)

_STRUCTURE_ALL = ("feature", "psm", "pg", "run", "sample")

# Above this consensusXML size, auto-select the streaming reader: the pyopenms
# in-memory load needs ~0.8x the file in RAM, so ~4 GB (≈3 GB RAM) is a safe point
# to switch to the low-memory path on a normal node.
_STREAM_THRESHOLD_BYTES = 4 * 1024**3


def _should_stream(path: str) -> bool:
    try:
        return Path(path).stat().st_size > _STREAM_THRESHOLD_BYTES
    except OSError:
        return False


def _validate_structures(structures: tuple[str, ...], sdrf_path: Optional[str]) -> None:
    """Reject unknown structure names and run/sample requested without an SDRF."""
    unknown = [s for s in structures if s not in _STRUCTURE_ALL]
    if unknown:
        raise ValueError(f"Unknown structure(s) {unknown}; valid values are {list(_STRUCTURE_ALL)}")
    metadata = set(structures).intersection({"run", "sample"})
    if metadata and not sdrf_path:
        raise ValueError(f"An SDRF is required to write {sorted(metadata)}")


def _write_sdrf_metadata(
    output_folder: Path,
    output_prefix: str,
    sdrf_path: str,
    requested: set[str],
) -> dict[str, Path]:
    """Write the requested SDRF-backed run/sample structures."""
    metadata = requested.intersection({"run", "sample"})
    if not metadata:
        return {}

    from qpx.converters.sdrf import SdrfConverter

    paths = {
        "sample": output_folder / f"{output_prefix}.sample.parquet",
        "run": output_folder / f"{output_prefix}.run.parquet",
    }
    with SdrfConverter() as sdrf_converter:
        sdrf_converter.convert(
            sdrf_path=sdrf_path,
            sample_output=str(paths["sample"]) if "sample" in metadata else None,
            run_output=str(paths["run"]) if "run" in metadata else None,
        )
    return {name: paths[name] for name in metadata}


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
        requested = set(structures)

        out = Path(output_folder)
        out.mkdir(parents=True, exist_ok=True)
        written: dict[str, Path] = {}

        # Parse the consensusXML once and share the map across the feature/psm/pg
        # adapters (they only read it). The in-memory pyopenms load dominates cost
        # and RAM; for large files use the streaming reader instead (re-parses per
        # adapter pass but never holds the whole map — output is identical).
        cm = None
        if {"feature", "psm", "pg"}.intersection(structures):
            use_stream = streaming if streaming is not None else _should_stream(consensusxml_path)
            if use_stream:
                from qpx.converters.openms_consensus.streaming import StreamingConsensusMap

                cm = StreamingConsensusMap(consensusxml_path)
            else:
                cm = load_consensus_map(consensusxml_path)
            # Ground-truth check: the channels read from the consensusXML maps must
            # match the SDRF comment[label] set — warn on any inconsistency (e.g. a
            # mis-declared plex or the wrong SDRF) rather than silently trusting one.
            if sdrf_path:
                for msg in check_channels_vs_sdrf(cm, sdrf_path):
                    _log.warning("consensusXML/SDRF channel mismatch: %s", msg)

        if "feature" in structures:
            # Share the protein-group leader map so feature.anchor_protein matches pg.
            anchor_map = accession_to_anchor(cm) if cm is not None else None
            recs = consensus_features_to_records(cm=cm, anchor_map=anchor_map)
            path = out / f"{output_prefix}.feature.parquet"
            with FeatureWriter(str(path), creator=creator) as w:
                if recs:
                    w.write_batch(recs)
            written["feature"] = path

        if "psm" in structures:
            recs = consensus_psms_to_records(cm=cm)
            path = out / f"{output_prefix}.psm.parquet"
            with PsmWriter(str(path), creator=creator) as w:
                if recs:
                    w.write_batch(recs)
            written["psm"] = path

        if "pg" in structures:
            recs = consensus_protein_groups_to_records(sdrf_path=sdrf_path, cm=cm, top=pg_top)
            path = out / f"{output_prefix}.pg.parquet"
            with PgWriter(str(path), creator=creator) as w:
                if recs:
                    w.write_batch(recs)
            written["pg"] = path

        written.update(_write_sdrf_metadata(out, output_prefix, sdrf_path, requested))

        return written
