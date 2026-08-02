"""QuantMS orchestrator — composes PSM, Feature, and PG adapters."""

from __future__ import annotations

import logging
from pathlib import Path

from qpx._version import __version__
from qpx.converters.mappings import get_tool_meta
from qpx.converters.mztab import (
    extract_modifications,
    load_msstats,
    load_mztab_sections,
)
from qpx.converters.orchestrator import BaseOrchestrator
from qpx.converters.quantms.feature_adapter import QuantmsFeatureAdapter
from qpx.converters.quantms.pg_adapter import QuantmsPgAdapter
from qpx.converters.quantms.psm_adapter import QuantmsPsmAdapter
from qpx.converters.sdrf import SdrfConverter
from qpx.core.constants import FEATURE, ONTOLOGY, PG, PSM, RUN, SAMPLE
from qpx.core.engine import create_converter_connection
from qpx.core.scores import (
    field_ontology_entries,
    modification_ontology_entries,
    score_ontology_entries,
)

logger = logging.getLogger(__name__)


class QuantMSConverter(BaseOrchestrator):
    """Orchestrate full QuantMS mzTab conversion to QPX format.

    Delegates to individual adapters based on the requested structures.
    """

    def __init__(
        self,
        mztab_path,
        sdrf_file,
        msstats_file=None,
        database_path=None,
        compression: str = "zstd",
    ):
        self.mztab_path = str(mztab_path)
        self.sdrf_file = str(sdrf_file)
        self.msstats_file = str(msstats_file) if msstats_file else None
        self.database_path = str(database_path) if database_path else None
        self._compression = compression
        self._resolved_mappings_by_view: dict[str, dict] = {}

    def convert(
        self,
        output_folder,
        output_prefix="quantms",
        structures=None,
        project_accession=None,
    ):
        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)

        if structures is None:
            structures = [PSM, FEATURE, PG]

        ontology_entries: list[dict] = []

        # Always convert SDRF for sample/run metadata
        with SdrfConverter() as sdrf_conv:
            sdrf_conv.convert(
                sdrf_path=self.sdrf_file,
                sample_output=str(output_folder / f"{output_prefix}.sample.parquet"),
                run_output=str(output_folder / f"{output_prefix}.run.parquet"),
            )
            ontology_entries.extend(sdrf_conv.run_ontology_entries())
            logger.info("SDRF conversion complete")

        # Load enzyme name from SDRF for missed-cleavage computation
        enzyme_name: str | None = None
        try:
            from qpx.core.sdrf import SDRFHandler

            handler = SDRFHandler(self.sdrf_file)
            enzymes = handler.get_enzymes()
            if enzymes:
                enzyme_name = str(enzymes[0])
        except Exception:
            logger.debug("Could not load enzyme from SDRF for missed cleavages")

        # Shared DuckDB connection — parse mzTab and MSstats once
        needs_mztab = any(s in structures for s in (PSM, FEATURE, PG))
        if not needs_mztab:
            self._write_dataset(
                output_folder,
                output_prefix,
                project_accession,
                software_name="OpenMS/quantms",
                software_version="quantms",
            )
            return

        shared_conn = create_converter_connection()
        provenance_records: list[dict] = []
        try:
            load_mztab_sections(shared_conn, self.mztab_path)
            if self.msstats_file:
                load_msstats(shared_conn, self.msstats_file)

            if PSM in structures:
                with QuantmsPsmAdapter(conn=shared_conn, compression=self._compression) as adapter:
                    adapter.convert(
                        mztab_path=self.mztab_path,
                        output_path=str(output_folder / f"{output_prefix}.psm.parquet"),
                        enzyme_name=enzyme_name,
                    )
                    ontology_entries.extend(score_ontology_entries(adapter.get_discovered_scores(), view=PSM))
                    self._resolved_mappings_by_view[PSM] = adapter.get_resolved_columns()
                    logger.info("PSM conversion complete")

            if FEATURE in structures and self.msstats_file:
                with QuantmsFeatureAdapter(conn=shared_conn, compression=self._compression) as adapter:
                    adapter.convert(
                        mztab_path=self.mztab_path,
                        msstats_path=self.msstats_file,
                        output_path=str(output_folder / f"{output_prefix}.feature.parquet"),
                        sdrf_path=self.sdrf_file,
                        enzyme_name=enzyme_name,
                    )
                    ontology_entries.extend(score_ontology_entries(adapter.get_discovered_scores(), view=FEATURE))
                    logger.info("Feature conversion complete")

            feature_path = output_folder / f"{output_prefix}.feature.parquet"
            if PG in structures and not feature_path.exists():
                raise ValueError(
                    "PG output was requested, but required feature input is missing: "
                    f"{feature_path}. Provide `msstats_file` and include `feature` "
                    "in structures (or pre-generate feature.parquet) before requesting `pg`."
                )

            if PG in structures:
                with QuantmsPgAdapter(conn=shared_conn, compression=self._compression) as adapter:
                    adapter.convert(
                        mztab_path=self.mztab_path,
                        feature_path=str(feature_path),
                        output_path=str(output_folder / f"{output_prefix}.pg.parquet"),
                        sdrf_path=self.sdrf_file,
                    )
                    ontology_entries.extend(score_ontology_entries(adapter.get_discovered_scores(), view=PG))
                    self._resolved_mappings_by_view[PG] = adapter.get_resolved_columns()
                    logger.info("PG conversion complete")

            # PTM ontology entries from mzTab modification metadata
            modifications_meta = extract_modifications(shared_conn)
            ontology_entries.extend(modification_ontology_entries(modifications_meta, view=PSM))

            # Field-level CV term entries with source provenance
            for view_name, mappings in self._resolved_mappings_by_view.items():
                ontology_entries.extend(
                    field_ontology_entries(
                        view=view_name,
                        resolved_mappings=mappings,
                        tool_name=get_tool_meta("quantms")["tool_name"],
                    )
                )

            # Build provenance from mzTab metadata
            provenance_records = self._build_provenance(shared_conn, self.mztab_path)
        finally:
            shared_conn.close()

        self._write_ontology(output_folder, output_prefix, ontology_entries)
        self._write_provenance(output_folder, output_prefix, provenance_records)
        self._write_dataset(
            output_folder,
            output_prefix,
            project_accession,
            software_name="OpenMS/quantms",
            software_version="quantms",
            provenance_records=provenance_records,
        )
        self._write_mudata(output_folder, output_prefix)

    # ------------------------------------------------------------------
    # Provenance extraction from mzTab metadata
    # ------------------------------------------------------------------

    @staticmethod
    def _build_provenance(
        conn,
        mztab_path: str | None = None,
    ) -> list[dict]:
        """Extract processing provenance from mzTab metadata.

        Reads ``software[N]`` entries from the metadata table to identify
        which search engines and tools were used.
        """
        import re

        provenance_records: list[dict] = []
        step_order = 0

        try:
            rows = conn.execute("""
                SELECT key, value FROM metadata
                WHERE key LIKE 'software[%'
                  AND key NOT LIKE '%-%'
                ORDER BY key
                """).fetchall()
        except Exception:
            rows = []

        for key, value in rows:
            step_order += 1
            # Parse CV term: [CV, accession, name, value]
            parts = value.replace("[", "").replace("]", "").split(",")
            tool_name = parts[2].strip() if len(parts) >= 3 else value.strip()

            # Try to extract version from the same metadata entry
            version = None
            if len(parts) >= 4 and parts[3].strip():
                version = parts[3].strip()

            # Determine step category from software index
            m = re.search(r"\[(\d+)\]", key)
            idx = int(m.group(1)) if m else step_order

            provenance_records.append(
                {
                    "step_order": idx,
                    "step_category": "database_search",
                    "step_name": "spectrum_identification",
                    "tool_name": tool_name,
                    "tool_version": version,
                    "tool_uri": None,
                    "parameters": None,
                    "config": None,
                    "output_views": [PSM, FEATURE, PG],
                }
            )

        # Add QPX conversion step
        step_order = max((r["step_order"] for r in provenance_records), default=0) + 1
        params = None
        if mztab_path:
            params = [{"key": "mztab_path", "value": Path(mztab_path).name}]
        provenance_records.append(
            {
                "step_order": step_order,
                "step_category": "format_conversion",
                "step_name": "mztab_to_qpx",
                "tool_name": "qpx",
                "tool_version": __version__,
                "tool_uri": None,
                "parameters": params,
                "config": None,
                "output_views": [PSM, FEATURE, PG, SAMPLE, RUN, ONTOLOGY],
            }
        )

        return provenance_records
