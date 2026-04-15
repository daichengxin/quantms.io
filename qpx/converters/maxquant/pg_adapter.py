"""MaxQuant PG adapter -- proteinGroups.txt to pg.parquet.

Loads ``proteinGroups.txt`` into DuckDB, transforms protein-group rows
into ``PgSchema``, and writes through ``PgWriter``.

Key schema changes:
    - ``reference_file_name`` -> ``run_file_name``
    - ``is_decoy`` (int32, 0/1) -> ``is_decoy`` (bool)
    - Intensities struct: ``{sample_accession, channel, intensity}`` ->
      ``{label, intensity}``
"""

from __future__ import annotations

import logging
import re
from typing import Optional

import pandas as pd

from qpx.converters.base import resolve_columns
from qpx.converters.mappings import get_field_mappings
from qpx.converters.maxquant.base_adapter import MaxQuantBaseAdapter
from qpx.converters.maxquant.constants import TMT_LABEL_TO_MQ_COL
from qpx.converters.utils import mq_flag_to_bool, safe_float
from qpx.writers.pg import PgWriter

logger = logging.getLogger(__name__)

# Derive field map from central mappings
_PG_MAP = get_field_mappings("maxquant", "pg")


class MaxQuantPgAdapter(MaxQuantBaseAdapter):
    """Convert MaxQuant ``proteinGroups.txt`` to ``pg.parquet``.

    Usage::

        with MaxQuantPgAdapter() as adapter:
            adapter.convert(
                protein_groups_path="proteinGroups.txt",
                output_path="pg.parquet",
                sdrf_path="sdrf.tsv",
            )
    """

    def convert(
        self,
        protein_groups_path: str,
        output_path: str,
        sdrf_path: Optional[str] = None,
        chunksize: int = 100_000,
        creator: str = "maxquant",
    ) -> None:
        """Run the proteinGroups.txt -> pg.parquet conversion.

        Args:
            protein_groups_path: Path to MaxQuant ``proteinGroups.txt``.
            output_path: Destination Parquet path.
            sdrf_path: Optional SDRF for sample mapping.
            chunksize: Rows per batch.
            creator: Creator tag in Parquet metadata.
        """
        # Step 1: Load proteinGroups.txt into DuckDB
        self._load_protein_groups(protein_groups_path)

        # Step 2: Resolve column mappings against actual input columns
        actual_cols = {
            c[0]
            for c in self._conn.execute(
                "SELECT column_name FROM information_schema.columns WHERE table_name='protein_groups'"
            ).fetchall()
        }
        self._resolved = resolve_columns(_PG_MAP, actual_cols)

        # Step 3: Load SDRF mapping
        _, experiment_type, tmt_channels = self._load_sdrf(sdrf_path)

        # Step 4: Detect intensity columns in the data
        intensity_cols = self._detect_intensity_columns(experiment_type)

        # Step 5: Stream and transform
        self.logger.info("Transforming MaxQuant protein groups ...")

        with PgWriter(output_path, creator=creator, compression=self._compression) as writer:
            for batch in self._query_batched("SELECT * FROM protein_groups", chunksize):
                df = batch.to_pandas()
                records = self._transform_batch(df, intensity_cols, tmt_channels)
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)

        self.logger.info(f"MaxQuant PG conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_protein_groups(self, path: str) -> None:
        """Load proteinGroups.txt into DuckDB."""
        self._conn.execute(
            "CREATE TABLE protein_groups AS SELECT * FROM read_csv_auto($1,"
            " delim='\\t', header=true, auto_detect=true, null_padding=true)",
            [path],
        )
        count = self._conn.execute("SELECT COUNT(*) FROM protein_groups").fetchone()[0]
        self.logger.info(f"Loaded {count:,} MaxQuant protein groups")

    def _detect_intensity_columns(self, experiment_type: str) -> dict[str, list[str]]:
        """Detect sample-specific intensity columns in proteinGroups.txt."""
        cols = self._conn.execute(
            "SELECT column_name FROM information_schema.columns WHERE table_name='protein_groups'"
        ).fetchall()
        col_names = [c[0] for c in cols]

        # Find Intensity <sample> columns (LFQ-style)
        intensity_cols = [c for c in col_names if c.startswith("Intensity ") and c != "Intensity"]
        lfq_cols = [c for c in col_names if c.startswith("LFQ intensity ")]
        ibaq_cols = [c for c in col_names if c.startswith("iBAQ ") and c != "iBAQ"]

        # TMT/iTRAQ: "Reporter intensity N <experiment>" columns
        exp_type_clean = re.sub(r"\d+", "", experiment_type).upper()
        reporter_cols: list[str] = []
        if exp_type_clean in ("TMT", "ITRAQ"):
            reporter_cols = [c for c in col_names if re.match(r"Reporter intensity \d+ .+", c) and "corrected" not in c.lower()]

        return {
            "intensity": intensity_cols,
            "lfq": lfq_cols,
            "ibaq": ibaq_cols,
            "reporter": reporter_cols,
        }

    # ------------------------------------------------------------------
    # Transform
    # ------------------------------------------------------------------

    def _transform_batch(
        self,
        df: pd.DataFrame,
        intensity_cols: dict[str, list[str]],
        tmt_channels: list[str] | None = None,
    ) -> list[dict]:
        records: list[dict] = []
        skipped = 0
        for row in df.to_dict("records"):
            try:
                recs = self._transform_row(row, intensity_cols, tmt_channels or [])
                records.extend(recs)
            except Exception as e:
                skipped += 1
                self.logger.debug(f"Skipping MaxQuant PG row: {e}")
        if skipped:
            total = skipped + len(records)
            self.logger.warning(
                "Skipped %d / %d rows (%.1f%%) in batch",
                skipped,
                total,
                100 * skipped / total if total else 0,
            )
        return records

    @staticmethod
    def _extract_gene_names(fasta_headers: str) -> list[str] | None:
        """Extract gene names from Fasta header lines using ``GN=`` tag."""
        genes: list[str] = []
        for header in fasta_headers.split(";"):
            m = re.search(r"GN=([^\s;]+)", header)
            if m:
                gene = m.group(1)
                if gene not in genes:
                    genes.append(gene)
        return genes or None

    def _transform_row(
        self,
        row,
        intensity_cols: dict[str, list[str]],
        tmt_channels: list[str] | None = None,
    ) -> list[dict]:
        """Transform a single proteinGroups.txt row.

        MaxQuant PG rows contain intensities for ALL runs in one row, so
        we expand them into one QPX PG record per run.
        """
        records: list[dict] = []
        r = self._resolved  # shorthand for resolved column mappings

        # Protein group identity
        pg_acc_cell = row.get(r.get("pg_accessions", "Protein IDs"))
        pg_acc_raw = str(pg_acc_cell) if pd.notna(pg_acc_cell) and pg_acc_cell is not None else ""
        pg_accessions = [self._norm_uniprot_id(a.strip()) for a in pg_acc_raw.split(";") if a.strip()] if pg_acc_raw else []

        pg_names_raw = row.get(r.get("pg_names", "Protein names"))
        pg_names = str(pg_names_raw).split(";") if pd.notna(pg_names_raw) and pg_names_raw else None
        # Fallback: use bare accessions as names (already normalised by _norm_uniprot_id)
        if not pg_names and pg_accessions:
            pg_names = list(pg_accessions)

        gg_raw = row.get(r.get("gg_accessions", "Gene names"))
        gg_accessions = str(gg_raw).split(";") if pd.notna(gg_raw) and gg_raw else None
        # Fallback: try Fasta headers for gene names
        if not gg_accessions:
            fasta_col = r.get("fasta_headers", "Fasta headers")
            fasta_raw = row.get(fasta_col)
            if pd.notna(fasta_raw) and fasta_raw:
                gg_accessions = self._extract_gene_names(str(fasta_raw))

        # Anchor protein (first in Majority protein IDs, or first PG accession)
        majority_raw = row.get(r.get("anchor_protein", "Majority protein IDs"))
        if pd.notna(majority_raw) and majority_raw:
            anchor_protein = self._norm_uniprot_id(str(majority_raw).split(";", maxsplit=1)[0])
        elif pg_accessions:
            anchor_protein = pg_accessions[0]
        else:
            anchor_protein = ""

        # Quality metrics
        global_qvalue = safe_float(row.get(r.get("global_qvalue", "Q-value")))
        is_decoy = mq_flag_to_bool(row.get(r.get("is_decoy", "Reverse"), ""))
        contaminant_val = mq_flag_to_bool(row.get(r.get("contaminant", "Potential contaminant"), ""))

        # Sequence coverage and molecular weight
        seq_coverage = safe_float(row.get(r.get("sequence_coverage", "Sequence coverage [%]")))
        mol_weight = safe_float(row.get(r.get("molecular_weight", "Mol. weight [kDa]")))

        # Peptide counts
        peptide_count_total = int(row.get(r.get("peptide_count_total", "Peptides"), 0) or 0)
        peptide_count_unique = int(row.get(r.get("peptide_count_unique", "Unique peptides"), 0) or 0)
        # Additional scores
        andromeda = safe_float(row.get(r.get("andromeda_score", "Score")))
        additional_scores = []
        if andromeda is not None:
            additional_scores.append(
                {
                    "score_name": "andromeda_score",
                    "score_value": andromeda,
                    "higher_better": True,
                }
            )

        # Peptides per protein
        peptides = [{"protein_name": acc, "peptide_count": peptide_count_total} for acc in pg_accessions]

        # Build shared record skeleton (fields independent of run/channel)
        def _make_rec(run_name: str, intensities: list, additional_intensities: list) -> dict:
            return {
                "pg_accessions": pg_accessions,
                "pg_names": pg_names,
                "gg_accessions": gg_accessions,
                "gg_names": gg_accessions,
                "gg_qvalue": None,
                "anchor_protein": anchor_protein,
                "run_file_name": run_name,
                "global_qvalue": global_qvalue,
                "pg_qvalue": global_qvalue,
                "intensities": intensities,
                "additional_intensities": additional_intensities or None,
                "is_decoy": is_decoy,
                "contaminant": contaminant_val,
                "peptides": peptides,
                "peptide_counts": {
                    "unique_sequences": peptide_count_unique,
                    "total_sequences": peptide_count_total,
                },
                "feature_counts": {
                    "unique_features": peptide_count_unique,
                    "total_features": peptide_count_total,
                },
                "sequence_coverage": seq_coverage,
                "molecular_weight": mol_weight,
                "additional_scores": additional_scores or None,
                "cv_params": None,
            }

        # ── TMT/iTRAQ: one record per experiment, channels as intensities list ──
        reporter_cols = intensity_cols.get("reporter", [])
        if reporter_cols:
            # Detect ri_offset from column names (0 vs 1 indexed)
            ri_offset = 0 if any(re.match(r"Reporter intensity 0 ", c) for c in reporter_cols) else 1
            # Group by experiment name
            experiments: dict[str, dict[int, str]] = {}  # exp -> {col_idx: col_name}
            for col in reporter_cols:
                m = re.match(r"Reporter intensity (\d+) (.+)", col)
                if m:
                    ch_num = int(m.group(1))
                    exp = m.group(2)
                    experiments.setdefault(exp, {})[ch_num] = col
            for exp, ch_cols in experiments.items():
                intensities: list[dict] = []
                additional_intensities: list[dict] = []
                channels = tmt_channels or []
                for seq_idx, channel_name in enumerate(channels):
                    col_idx = TMT_LABEL_TO_MQ_COL.get(str(channel_name).strip().upper())
                    if col_idx is None:
                        col_idx = seq_idx
                    mq_col_num = col_idx + ri_offset
                    col_name = ch_cols.get(mq_col_num)
                    if not col_name:
                        continue
                    val = safe_float(row.get(col_name))
                    if val and val > 0:
                        intensities.append({"label": channel_name, "intensity": float(val)})
                        corr_col = col_name.replace("Reporter intensity", "Reporter intensity corrected")
                        corr_val = safe_float(row.get(corr_col))
                        if corr_val is not None:
                            additional_intensities.append(
                                {
                                    "label": channel_name,
                                    "intensities": [
                                        {"intensity_name": "corrected_reporter_intensity", "intensity_value": float(corr_val)}
                                    ],
                                }
                            )
                if intensities:
                    records.append(_make_rec(exp, intensities, additional_intensities))
        else:
            # ── LFQ: one record per run using "Intensity <run>" columns ──
            for intensity_col in intensity_cols.get("intensity", []):
                run_name = intensity_col.removeprefix("Intensity ")
                intensity_val = safe_float(row.get(intensity_col)) or 0.0
                if intensity_val <= 0:
                    continue
                lfq_val = safe_float(row.get(f"LFQ intensity {run_name}"))
                ibaq_val = safe_float(row.get(f"iBAQ {run_name}"))
                extra_vals = []
                if lfq_val is not None:
                    extra_vals.append({"intensity_name": "lfq", "intensity_value": float(lfq_val)})
                if ibaq_val is not None:
                    extra_vals.append({"intensity_name": "ibaq", "intensity_value": float(ibaq_val)})
                add_int = [{"label": "LFQ", "intensities": extra_vals}] if extra_vals else []
                records.append(
                    _make_rec(
                        run_name,
                        [{"label": "LFQ", "intensity": float(intensity_val)}],
                        add_int,
                    )
                )

        # If no per-sample intensity columns, emit one record with total Intensity
        if not records:
            total_intensity = safe_float(row.get(r.get("intensity", "Intensity"))) or 0.0
            if total_intensity > 0:
                records.append(
                    {
                        "pg_accessions": pg_accessions,
                        "pg_names": pg_names,
                        "gg_accessions": gg_accessions,
                        "gg_names": gg_accessions,  # Gene symbols serve as both accession and name
                        "gg_qvalue": None,
                        "anchor_protein": anchor_protein,
                        "run_file_name": "unknown",
                        "global_qvalue": global_qvalue,
                        "pg_qvalue": global_qvalue,
                        "intensities": [{"label": "LFQ", "intensity": float(total_intensity)}],
                        "additional_intensities": None,
                        "is_decoy": is_decoy,
                        "contaminant": contaminant_val,
                        "peptides": peptides,
                        "peptide_counts": {
                            "unique_sequences": peptide_count_unique,
                            "total_sequences": peptide_count_total,
                        },
                        "feature_counts": {
                            "unique_features": peptide_count_unique,
                            "total_features": peptide_count_total,
                        },
                        "sequence_coverage": seq_coverage,
                        "molecular_weight": mol_weight,
                        "additional_scores": additional_scores or None,
                        "cv_params": None,
                    }
                )

        return records
