"""Orchestrate consensusXML (+ SDRF) -> a QPX dataset (feature/psm/pg [+ run/sample]).

Interim quantms path while OpenMS ``-out_qpx`` is pre-1.1. pg is identification-
only (no protein intensity). run/sample come from the SDRF (reusing
:class:`SdrfConverter`) when an SDRF is provided.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from qpx.converters.openms_consensus.feature_adapter import consensus_features_to_records, load_consensus_map
from qpx.converters.openms_consensus.pg_adapter import consensus_protein_groups_to_records
from qpx.converters.openms_consensus.psm_adapter import consensus_psms_to_records
from qpx.writers.feature import FeatureWriter
from qpx.writers.pg import PgWriter
from qpx.writers.psm import PsmWriter

_STRUCTURE_ALL = ("feature", "psm", "pg", "run", "sample")


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
    ) -> dict[str, Path]:
        """Write the requested QPX views and return ``{structure: parquet path}``.

        ``structures`` selects which of feature/psm/pg/run/sample to emit; pg is
        identification-only (null protein intensity) in this interim path.
        """
        requested = set(structures)
        unknown = requested.difference(_STRUCTURE_ALL)
        if unknown:
            raise ValueError(f"Unknown OpenMS consensus structure(s): {sorted(unknown)}")
        metadata = requested.intersection({"run", "sample"})
        if metadata and not sdrf_path:
            raise ValueError(f"An SDRF is required to write {sorted(metadata)}")

        out = Path(output_folder)
        out.mkdir(parents=True, exist_ok=True)
        written: dict[str, Path] = {}

        # Parse the consensusXML once and share the in-memory map across the
        # feature/psm/pg adapters (they only read it) — the load dominates the
        # convert cost, so re-parsing per adapter would triple it.
        cm = None
        if {"feature", "psm", "pg"}.intersection(structures):
            cm = load_consensus_map(consensusxml_path)

        if "feature" in structures:
            recs = consensus_features_to_records(cm=cm)
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
            recs = consensus_protein_groups_to_records(sdrf_path=sdrf_path, cm=cm)
            path = out / f"{output_prefix}.pg.parquet"
            with PgWriter(str(path), creator=creator) as w:
                if recs:
                    w.write_batch(recs)
            written["pg"] = path

        written.update(_write_sdrf_metadata(out, output_prefix, sdrf_path, requested))

        return written
