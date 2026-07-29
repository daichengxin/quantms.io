"""CDAP orchestrator -- composes PSM, Feature, and PG adapters."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable, Optional

from qpx._version import __version__
from qpx.converters.cdap.feature_adapter import CdapFeatureAdapter
from qpx.converters.cdap.pg_adapter import CdapPgAdapter
from qpx.converters.cdap.psm_adapter import CdapPsmAdapter
from qpx.converters.mappings import get_tool_meta
from qpx.converters.orchestrator import BaseOrchestrator
from qpx.core.constants import FEATURE, ONTOLOGY, PG, PSM
from qpx.core.engine import create_converter_connection
from qpx.core.scores import field_ontology_entries, score_ontology_entries

logger = logging.getLogger(__name__)


class CdapConverter(BaseOrchestrator):
    """Orchestrate full CDAP conversion to QPX format.

    A CPTAC study directory (e.g. ``PDC000440``) contains many ``*.psm``
    tab-separated files.  The orchestrator runs three adapters in sequence:

    1. :class:`CdapPsmAdapter` -- one PSM record per ``.psm`` row.
    2. :class:`CdapFeatureAdapter` -- aggregate PSMs into peptidoform features.
    3. :class:`CdapPgAdapter` -- aggregate features into protein groups,
       reading the feature parquet just written.
    """

    def __init__(
        self,
        max_memory: str = "16GB",
        max_cpus: int = 4,
        compression: str = "zstd",
    ):
        self._max_memory = max_memory
        self._max_cpus = max_cpus
        self._compression = compression
        self._resolved_mappings_by_view: dict[str, dict] = {}

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @staticmethod
    def supported_structures() -> tuple[str, ...]:
        """Return the QPX structures this converter can produce."""
        return ("psm", "feature", "pg")

    def convert(
        self,
        psm_dir: str | Path,
        output_folder: str | Path,
        output_prefix: str = "cdap",
        structures: Optional[Iterable[str]] = None,
        batch_size: int = 200_000,
        project_accession: Optional[str] = None,
    ) -> None:
        """Convert a single CDAP study directory to QPX Parquet files.

        Args:
            psm_dir: Directory containing one study's ``*.psm`` files.
            output_folder: Destination directory for QPX outputs.
            output_prefix: Prefix for output filenames (default: ``cdap``).
            structures: Iterable of structures to produce (subset of
                ``{psm, feature, pg}``). When ``None`` all three are produced.
            batch_size: Rows per DuckDB->pandas batch.
            project_accession: PDC study accession written into
                ``dataset.parquet``.
        """
        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)
        prefix = output_prefix or "cdap"

        requested = self._normalize_structures(structures)
        logger.info("CDAP converter producing: %s", sorted(requested))

        ontology_entries: list[dict] = []
        psm_path = output_folder / f"{prefix}.psm.parquet"
        feature_path = output_folder / f"{prefix}.feature.parquet"
        pg_path = output_folder / f"{prefix}.pg.parquet"

        self._convert_psm_and_feature(requested, psm_dir, psm_path, feature_path, batch_size, ontology_entries)
        self._convert_pg(requested, feature_path, pg_path, batch_size, ontology_entries)
        self._write_metadata(output_folder, prefix, requested, ontology_entries, project_accession)

    # ------------------------------------------------------------------
    # Conversion sub-steps
    # ------------------------------------------------------------------

    def _convert_psm_and_feature(self, requested, psm_dir, psm_path, feature_path, batch_size, ontology_entries):
        """Run PSM and Feature adapters with shared DuckDB connection."""
        need_shared = PSM in requested and FEATURE in requested
        shared_conn = create_converter_connection(memory_limit=self._max_memory, threads=self._max_cpus) if need_shared else None
        try:
            if PSM in requested:
                with CdapPsmAdapter(
                    duckdb_memory=self._max_memory,
                    duckdb_threads=self._max_cpus,
                    compression=self._compression,
                    conn=shared_conn,
                ) as adapter:
                    adapter.convert(psm_dir=str(psm_dir), output_path=str(psm_path), chunksize=batch_size)
                    ontology_entries.extend(score_ontology_entries(adapter.get_discovered_scores(), view=PSM))
                    self._resolved_mappings_by_view[PSM] = adapter.get_resolved_columns()
                logger.info("CDAP PSM conversion complete")

            if FEATURE in requested:
                with CdapFeatureAdapter(
                    duckdb_memory=self._max_memory,
                    duckdb_threads=self._max_cpus,
                    compression=self._compression,
                    conn=shared_conn,
                ) as adapter:
                    adapter.convert(psm_dir=str(psm_dir), output_path=str(feature_path), chunksize=batch_size)
                    ontology_entries.extend(score_ontology_entries(adapter.get_discovered_scores(), view=FEATURE))
                    self._resolved_mappings_by_view[FEATURE] = adapter.get_resolved_columns()
                logger.info("CDAP feature conversion complete")
        finally:
            if shared_conn is not None:
                shared_conn.close()

    def _convert_pg(self, requested, feature_path, pg_path, batch_size, ontology_entries):
        """Run PG adapter using existing feature.parquet."""
        if PG not in requested:
            return
        if not feature_path.exists():
            logger.warning(
                "PG conversion requires feature.parquet but %s is missing; "
                "either request 'feature' as well or supply an existing file.",
                feature_path,
            )
            return
        with CdapPgAdapter(
            duckdb_memory=self._max_memory,
            duckdb_threads=self._max_cpus,
            compression=self._compression,
        ) as adapter:
            adapter.convert(feature_path=str(feature_path), output_path=str(pg_path), chunksize=batch_size)
            ontology_entries.extend(score_ontology_entries(adapter.get_discovered_scores(), view=PG))
        logger.info("CDAP PG conversion complete")

    def _write_metadata(self, output_folder, prefix, requested, ontology_entries, project_accession):
        """Write ontology, provenance, and dataset files."""
        tool_meta = get_tool_meta("cdap")
        for view_name, mappings in self._resolved_mappings_by_view.items():
            ontology_entries.extend(
                field_ontology_entries(view=view_name, resolved_mappings=mappings, tool_name=tool_meta["tool_name"])
            )
        self._write_ontology(output_folder, prefix, ontology_entries)
        self._write_provenance(output_folder, prefix, sorted(requested))
        self._write_dataset(
            output_folder,
            prefix,
            project_accession,
            software_name=tool_meta["tool_name"],
            software_version=tool_meta.get("tool_versions"),
        )

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _normalize_structures(structures: Optional[Iterable[str]]) -> set[str]:
        """Resolve the requested structures (defaulting to PSM+Feature+PG)."""
        valid = {PSM, FEATURE, PG}
        if structures is None:
            return set(valid)
        requested = {s.strip().lower() for s in structures if s and s.strip()}
        unknown = requested - valid
        if unknown:
            raise ValueError(f"Unknown CDAP structures: {sorted(unknown)}")
        return requested or set(valid)

    def _write_provenance(
        self,
        output_folder: Path,
        prefix: str,
        records: list[str],
    ) -> None:
        """Write provenance.parquet with CDAP + QPX conversion steps."""
        prov_records = [
            {
                "step_order": 1,
                "step_category": "database_search",
                "step_name": "cdap_msgf_search",
                "tool_name": "MS-GF+",
                "tool_version": None,
                "tool_uri": "https://github.com/MSGFPlus/msgfplus",
                "parameters": None,
                "config": None,
                "output_views": records,
            },
            {
                "step_order": 2,
                "step_category": "format_conversion",
                "step_name": "cdap_to_qpx",
                "tool_name": "qpx",
                "tool_version": __version__,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                # CDAP itself produces psm/feature/pg + ontology; sample/run come
                # from PDC metadata in pdc2qpx, not from this converter.
                "output_views": records + [ONTOLOGY],
            },
        ]
        super()._write_provenance(output_folder, prefix, prov_records)
