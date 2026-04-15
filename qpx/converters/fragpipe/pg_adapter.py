"""FragPipe PG adapter -- combined_protein.tsv to pg.parquet.

Loads a FragPipe ``combined_protein.tsv`` into DuckDB, transforms
protein-group rows into ``PgSchema``, and writes through ``PgWriter``.

Key schema changes:
    - ``is_decoy`` (int, 0/1) -> ``is_decoy`` (bool)
    - Intensities struct: ``{sample_accession, channel, intensity}`` ->
      ``{label, intensity}``
"""

from __future__ import annotations

import logging

import pandas as pd

from qpx.converters.base import BaseConverter, resolve_columns
from qpx.converters.mappings import get_field_mappings
from qpx.converters.utils import safe_float
from qpx.core.sql import escape_path, sql_build
from qpx.writers.pg import PgWriter

logger = logging.getLogger(__name__)

# Derive field map from central YAML mappings
_PG_MAP = get_field_mappings("fragpipe", "pg")


class FragPipePgAdapter(BaseConverter):
    """Convert FragPipe ``combined_protein.tsv`` to ``pg.parquet``.

    Usage::

        with FragPipePgAdapter() as adapter:
            adapter.convert(
                protein_path="combined_protein.tsv",
                output_path="pg.parquet",
            )
    """

    def convert(
        self,
        protein_path: str,
        output_path: str,
        chunksize: int = 100_000,
        creator: str = "fragpipe",
    ) -> None:
        """Run the combined_protein.tsv -> pg.parquet conversion.

        Args:
            protein_path: Path to FragPipe ``combined_protein.tsv``.
            output_path: Destination Parquet path.
            chunksize: Rows per batch.
            creator: Creator tag in Parquet metadata.
        """
        # Step 1: Load protein file into DuckDB
        self._load_protein_file(protein_path)

        # Step 2: Resolve column mappings against actual input columns
        actual_cols = {
            c[0]
            for c in self._conn.execute(
                "SELECT column_name FROM information_schema.columns WHERE table_name='fragpipe_proteins'"
            ).fetchall()
        }
        self._resolved = resolve_columns(_PG_MAP, actual_cols)

        # Step 3: Detect per-experiment intensity columns
        experiment_cols = self._detect_experiment_columns()

        # Step 3: Stream and transform
        self.logger.info("Transforming FragPipe protein groups ...")

        with PgWriter(output_path, creator=creator, compression=self._compression) as writer:
            for batch in self._query_batched("SELECT * FROM fragpipe_proteins", chunksize):
                df = batch.to_pandas()
                records = self._transform_batch(df, experiment_cols)
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)

        self.logger.info(f"FragPipe PG conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_protein_file(self, path: str) -> None:
        """Load combined_protein.tsv into DuckDB."""
        self._conn.execute(
            sql_build(
                """CREATE TABLE fragpipe_proteins AS
            SELECT * FROM read_csv_auto('$path',
                delim='\t', header=true, auto_detect=true)""",
                path=escape_path(path),
            )
        )
        count = self._conn.execute("SELECT COUNT(*) FROM fragpipe_proteins").fetchone()[0]
        self.logger.info(f"Loaded {count:,} FragPipe protein group rows")

    def _detect_experiment_columns(self) -> list[str]:
        """Detect per-experiment intensity columns.

        FragPipe ``combined_protein.tsv`` has columns like:
            ``<experiment> Total Intensity``
            ``<experiment> Spectral Count``
            ``<experiment> Unique Spectral Count``
        """
        cols = self._conn.execute(
            "SELECT column_name FROM information_schema.columns WHERE table_name='fragpipe_proteins'"
        ).fetchall()
        col_names = [c[0] for c in cols]

        # Find columns ending with " Total Intensity" or " MaxLFQ Intensity"
        intensity_cols = [c for c in col_names if c.endswith(" Total Intensity") or c.endswith(" MaxLFQ Intensity")]

        # Extract experiment names
        experiments = set()
        for col in intensity_cols:
            for suffix in [" Total Intensity", " MaxLFQ Intensity"]:
                if col.endswith(suffix):
                    experiments.add(col[: -len(suffix)])

        return sorted(experiments)

    # ------------------------------------------------------------------
    # Transform
    # ------------------------------------------------------------------

    def _transform_batch(self, df: pd.DataFrame, experiments: list[str]) -> list[dict]:
        records: list[dict] = []
        skipped = 0
        for row in df.to_dict("records"):
            try:
                recs = self._transform_row(row, experiments)
                records.extend(recs)
            except Exception as e:
                skipped += 1
                self.logger.debug(f"Skipping FragPipe PG row: {e}")
        if skipped:
            total = skipped + len(records)
            self.logger.warning(
                "Skipped %d / %d rows (%.1f%%) in batch",
                skipped,
                total,
                100 * skipped / total if total else 0,
            )
        return records

    def _transform_row(self, row, experiments: list[str]) -> list[dict]:
        """Transform a single combined_protein.tsv row.

        Expands into one record per experiment (run).
        """
        records: list[dict] = []
        r = self._resolved  # shorthand for resolved column mappings

        # Protein group identity
        protein_raw = str(row.get(r.get("pg_accessions", "Protein"), ""))
        pg_accessions = protein_raw.split(",") if protein_raw else []

        gene_raw = row.get(r.get("gg_accessions", "Gene"), "")
        gg_accessions = str(gene_raw).split(",") if pd.notna(gene_raw) and gene_raw else None

        # Protein description/name
        description = row.get(r.get("pg_names", "Description"), "")
        pg_names = None
        if pd.notna(description) and description:
            pg_names = [str(description)]

        anchor_protein = pg_accessions[0] if pg_accessions else ""

        # Is decoy: detected from rev_/DECOY_/decoy_ prefix on any accession
        is_decoy = any(acc.strip().startswith(("rev_", "DECOY_", "decoy_")) for acc in pg_accessions) if pg_accessions else False

        # Combined counts
        total_peptides = int(row.get(r.get("peptide_count_total", "Combined Total Peptides"), 0) or 0)
        unique_peptides = int(row.get(r.get("peptide_count_unique", "Combined Unique Peptides"), 0) or 0)

        # Sequence coverage
        seq_coverage = safe_float(row.get(r.get("sequence_coverage", "Percent Coverage")))

        # Molecular weight
        mol_weight = safe_float(row.get(r.get("molecular_weight", "Protein Molecular Weight (Da)")))
        if mol_weight is not None:
            mol_weight = mol_weight / 1000.0  # Convert Da to kDa

        # Protein group q-value (Protein Probability, Protein FDR, etc.)
        pg_qvalue = safe_float(row.get(r.get("pg_qvalue", "Protein Probability")))

        # Peptides per protein
        peptides = [{"protein_name": acc, "peptide_count": total_peptides} for acc in pg_accessions]

        # Expand per experiment
        for experiment in experiments:
            total_int_col = f"{experiment} Total Intensity"
            maxlfq_col = f"{experiment} MaxLFQ Intensity"
            spec_count_col = f"{experiment} Spectral Count"
            unique_spec_col = f"{experiment} Unique Spectral Count"

            total_intensity = safe_float(row.get(total_int_col)) or 0.0
            if total_intensity <= 0:
                continue

            # Intensities (new schema: {label, intensity})
            label = "LFQ"
            intensities = [{"label": label, "intensity": float(total_intensity)}]

            # Additional intensities pre-computed by FragPipe (MaxLFQ)
            additional_intensities = []
            extra_vals = []
            maxlfq_val = safe_float(row.get(maxlfq_col))
            if maxlfq_val is not None:
                extra_vals.append({"intensity_name": "maxlfq", "intensity_value": float(maxlfq_val)})
            if extra_vals:
                additional_intensities.append({"label": label, "intensities": extra_vals})

            # Additional scores
            additional_scores = []
            spec_val = safe_float(row.get(spec_count_col))
            if spec_val is not None:
                additional_scores.append(
                    {
                        "score_name": "spectral_count",
                        "score_value": spec_val,
                        "higher_better": True,
                    }
                )

            # Per-experiment peptide/feature counts
            exp_unique = safe_float(row.get(unique_spec_col)) or 0
            exp_total = safe_float(row.get(spec_count_col)) or 0

            rec = {
                "pg_accessions": pg_accessions,
                "pg_names": pg_names,
                "gg_accessions": gg_accessions,
                "gg_names": gg_accessions,  # Gene symbols serve as both accession and name
                "gg_qvalue": None,
                "anchor_protein": anchor_protein,
                "run_file_name": experiment,
                "global_qvalue": None,
                "pg_qvalue": pg_qvalue,
                "intensities": intensities,
                "additional_intensities": additional_intensities or None,
                "is_decoy": is_decoy,
                "contaminant": None,
                "peptides": peptides,
                "peptide_counts": {
                    "unique_sequences": unique_peptides,
                    "total_sequences": total_peptides,
                },
                "feature_counts": {
                    "unique_features": int(exp_unique),
                    "total_features": int(exp_total),
                },
                "sequence_coverage": seq_coverage,
                "molecular_weight": mol_weight,
                "additional_scores": additional_scores or None,
                "cv_params": None,
            }
            records.append(rec)

        return records
