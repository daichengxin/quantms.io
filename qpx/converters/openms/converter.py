"""OpenMS QPX converter — enriches ProteomicsLFQ ``-out_qpx`` output.

OpenMS ``ProteomicsLFQ -out_qpx`` already produces QPX-compliant
``psm.parquet``, ``feature.parquet``, and ``pg.parquet``.  This converter
validates those files, copies them to the output folder, and generates the
missing metadata tables (``run``, ``sample``, ``ontology``, ``provenance``,
``dataset``) from an SDRF file.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from qpx._version import __version__
from qpx.converters.orchestrator import BaseOrchestrator
from qpx.converters.sdrf import SdrfConverter
from qpx.core.constants import FEATURE, ONTOLOGY, PG, PSM, RUN, SAMPLE
from qpx.core.data.loader import load_schema
from qpx.core.scores import score_ontology_entries

logger = logging.getLogger(__name__)

# Mapping from QPX view name to the schema loader name
_VIEW_SCHEMAS = {
    PSM: "psm",
    FEATURE: "feature",
    PG: "pg",
}


def _discover_parquet(qpx_dir: Path) -> dict[str, Path]:
    """Discover QPX parquet files in a directory.

    Looks for files matching ``*.psm.parquet``, ``*.feature.parquet``,
    and ``*.pg.parquet``.  Returns a dict mapping view name to path.
    """
    found: dict[str, Path] = {}
    for view in (PSM, FEATURE, PG):
        matches = sorted(qpx_dir.glob(f"*.{view}.parquet"))
        if matches:
            found[view] = matches[0]
            if len(matches) > 1:
                logger.warning(
                    "Multiple %s parquet files found in %s; using %s",
                    view,
                    qpx_dir,
                    matches[0].name,
                )
    return found


def _collect_score_names(table_path: Path) -> set[str]:
    """Read ``additional_scores`` from a parquet file and collect score names."""
    table = pq.read_table(str(table_path), columns=["additional_scores"])
    names: set[str] = set()
    col = table.column("additional_scores")
    for row in col.to_pylist():
        if row:
            for score in row:
                if score and "name" in score and score["name"]:
                    names.add(score["name"])
    return names


def _validate_core(discovered: dict[str, Path]) -> None:
    """Validate each discovered parquet file against its QPX schema."""
    for view, path in discovered.items():
        schema = load_schema(_VIEW_SCHEMAS[view])
        table = pq.read_table(str(path))
        result = schema.validate_full(table)
        if not result.is_valid:
            errors = "; ".join(i.message for i in result.errors)
            raise ValueError(f"Validation failed for {path.name}: {errors}")
        logger.info(
            "%s: %s (%d rows, %d errors, %d warnings)",
            view,
            path.name,
            table.num_rows,
            len(result.errors),
            len(result.warnings),
        )


def _copy_core(
    discovered: dict[str, Path],
    output_folder: Path,
    output_prefix: str,
) -> dict[str, Path]:
    """Copy core parquet files to output directory, skipping same-location."""
    output_paths: dict[str, Path] = {}
    for view, src_path in discovered.items():
        dst = output_folder / f"{output_prefix}.{view}.parquet"
        if src_path.resolve() != dst.resolve():
            shutil.copy2(str(src_path), str(dst))
            logger.info("Copied %s -> %s", src_path.name, dst.name)
        else:
            logger.info("Skipping copy for %s (same location)", view)
        output_paths[view] = dst
    return output_paths


class OpenMSConverter(BaseOrchestrator):
    """Orchestrate OpenMS QPX enrichment to a full QPX dataset.

    Usage::

        converter = OpenMSConverter(
            qpx_dir="path/to/openms_qpx_output",
            sdrf_path="metadata.sdrf.tsv",
        )
        converter.convert(
            output_folder="./qpx_full",
            output_prefix="openms",
        )
    """

    def __init__(
        self,
        qpx_dir: str | Path,
        sdrf_path: str | Path | None = None,
        compression: str = "zstd",
    ):
        """Initialize the OpenMS QPX converter.

        Parameters
        ----------
        qpx_dir : str or Path
            Directory containing OpenMS ``-out_qpx`` parquet files.
        sdrf_path : str, Path or None
            Optional SDRF metadata file for sample/run generation.
        compression : str
            Parquet compression codec (default ``zstd``).

        """
        self.qpx_dir = Path(qpx_dir)
        self.sdrf_path = str(sdrf_path) if sdrf_path else None
        self._compression = compression

    def convert(
        self,
        output_folder: str | Path,
        output_prefix: str = "openms",
        project_accession: str | None = None,
    ) -> None:
        """Run the full enrichment pipeline.

        1. Discover and validate core QPX parquet files
        2. Copy core files to output folder (if needed)
        3. Convert SDRF to run.parquet + sample.parquet
        4. Collect score names for ontology
        5. Write ontology.parquet, provenance.parquet, dataset.parquet
        """
        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)

        discovered = self.discover_and_validate()
        output_paths = _copy_core(discovered, output_folder, output_prefix)

        ontology_entries = self._convert_sdrf(output_folder, output_prefix)
        ontology_entries.extend(self._collect_ontology(output_paths))

        self._write_ontology(output_folder, output_prefix, ontology_entries)
        structures = list(discovered.keys())
        provenance_records = self._build_provenance(structures)
        self._write_provenance(output_folder, output_prefix, provenance_records)
        self._write_dataset(
            output_folder,
            output_prefix,
            project_accession,
            software_name="OpenMS/ProteomicsLFQ",
            software_version=None,
            provenance_records=provenance_records,
        )
        logger.info("OpenMS QPX enrichment complete -> %s", output_folder)

    def discover_and_validate(self) -> dict[str, Path]:
        """Discover and validate core QPX parquet files."""
        discovered = _discover_parquet(self.qpx_dir)
        if not discovered:
            raise FileNotFoundError(
                f"No QPX parquet files (*.psm.parquet, *.feature.parquet, *.pg.parquet) found in {self.qpx_dir}"
            )
        logger.info(
            "Discovered %d core QPX file(s): %s",
            len(discovered),
            ", ".join(f"{v}={p.name}" for v, p in discovered.items()),
        )
        _validate_core(discovered)
        return discovered

    def _convert_sdrf(
        self,
        output_folder: Path,
        output_prefix: str,
    ) -> list[dict]:
        """Convert SDRF to run + sample, returning ontology entries."""
        if not self.sdrf_path:
            logger.warning("No SDRF path provided — skipping sample/run conversion")
            return []
        with SdrfConverter(compression=self._compression) as sdrf_conv:
            sdrf_conv.convert(
                sdrf_path=self.sdrf_path,
                sample_output=str(output_folder / f"{output_prefix}.sample.parquet"),
                run_output=str(output_folder / f"{output_prefix}.run.parquet"),
            )
            entries = list(sdrf_conv.run_ontology_entries())
        logger.info("SDRF conversion complete (sample + run)")
        return entries

    @staticmethod
    def _collect_ontology(output_paths: dict[str, Path]) -> list[dict]:
        """Collect score-based ontology entries from core tables."""
        entries: list[dict] = []
        all_scores: set[str] = set()
        for view, path in output_paths.items():
            try:
                scores = _collect_score_names(path)
            except (KeyError, pa.ArrowInvalid):
                logger.debug("No additional_scores in %s", path.name)
                continue
            if scores:
                entries.extend(score_ontology_entries(scores, view=view))
                all_scores |= scores
        if all_scores:
            logger.info("Discovered scores: %s", sorted(all_scores))
        return entries

    @staticmethod
    def _build_provenance(structures: list[str]) -> list[dict]:
        """Build provenance records for OpenMS + QPX conversion."""
        return [
            {
                "step_order": 1,
                "step_category": "quantification",
                "step_name": "label_free_quantification",
                "tool_name": "OpenMS/ProteomicsLFQ",
                "tool_version": None,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                "output_views": structures,
            },
            {
                "step_order": 2,
                "step_category": "format_conversion",
                "step_name": "openms_qpx_enrichment",
                "tool_name": "qpx",
                "tool_version": __version__,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                "output_views": [SAMPLE, RUN, ONTOLOGY],
            },
        ]
