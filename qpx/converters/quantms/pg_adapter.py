"""QuantMS PG adapter -- aggregate feature.parquet to pg.parquet.

Reads a feature.parquet file (produced by QuantmsFeatureAdapter) and the
mzTab proteins table, aggregates feature-level intensities per
(anchor_protein, run_file_name) group, and writes the results through
``PgWriter``.
"""

from __future__ import annotations

import logging
import re
from typing import Optional

import numpy as np
import pandas as pd

from qpx.converters.base import BaseConverter
from qpx.converters.mztab import load_mztab_sections
from qpx.converters.utils import safe_float
from qpx.core.sql import escape_path, sql_build
from qpx.writers.pg import PgWriter

logger = logging.getLogger(__name__)


class QuantmsPgAdapter(BaseConverter):
    """Convert QuantMS features + mzTab proteins to ``pg.parquet``.

    Reads the already-produced feature.parquet and the mzTab proteins
    table to create protein-group level quantification.

    Usage::

        with QuantmsPgAdapter() as adapter:
            adapter.convert(
                mztab_path="results.mzTab",
                feature_path="output.feature.parquet",
                output_path="output.pg.parquet",
            )
    """

    def convert(
        self,
        mztab_path: str,
        feature_path: str,
        output_path: str,
        creator: str = "quantms",
        query_batch_size: int = 50_000,
        write_batch_size: int = 10_000,
        topn: Optional[int] = None,
        aggregation: str = "sum",
    ) -> None:
        """Run the feature.parquet + mzTab proteins -> pg.parquet conversion.

        Args:
            mztab_path: Path to the mzTab file (for protein group metadata).
            feature_path: Path to the feature.parquet file.
            output_path: Destination Parquet path.
            creator: Creator tag in Parquet metadata.
            topn: Number of top peptides to use for additional intensity
                calculation.  ``None`` (default) uses all peptides.
            aggregation: Aggregation method for additional intensities:
                ``"sum"`` (default) or ``"mean"``.
        """
        self._topn = topn
        self._aggregation = aggregation

        # Step 1: Load mzTab protein group metadata (two-pass: single + group)
        if not self._table_exists("proteins"):
            load_mztab_sections(self._conn, mztab_path)
        single_meta, group_meta = self._build_protein_meta()

        # Step 2: Stream sorted feature rows via DuckDB to avoid full materialization.
        safe_fp = escape_path(feature_path)
        self.logger.info(f"Reading features from {feature_path}")
        total_features = self._conn.execute(
            sql_build(
                "SELECT COUNT(*) FROM read_parquet('$path')",
                path=safe_fp,
            )
        ).fetchone()[0]
        self.logger.info(f"Loaded {total_features} features")

        # Step 3: Aggregate features per (anchor_protein, run_file_name) while streaming.
        self.logger.info("Aggregating protein groups from features ...")
        grouped_sql = sql_build(
            "SELECT * FROM read_parquet('$path') ORDER BY anchor_protein, run_file_name",
            path=safe_fp,
        )

        records_buffer: list[dict] = []
        total_records = 0
        current_anchor: Optional[str] = None
        current_run: Optional[str] = None
        current_rows: list[dict] = []
        total_groups = 0
        processed_groups = 0
        skipped_groups = 0
        failed_groups = 0
        null_key_count = 0
        max_bad_group_ratio = 0.20
        min_groups_for_ratio_failure = 5

        def _flush_current_group() -> None:
            nonlocal current_rows, total_records
            nonlocal total_groups, processed_groups, skipped_groups, failed_groups
            if not current_rows or current_anchor is None or current_run is None:
                return
            total_groups += 1
            try:
                rec = self._build_single_pg(
                    anchor_protein=current_anchor,
                    run_file_name=current_run,
                    features=pd.DataFrame.from_records(current_rows),
                    single_meta=single_meta,
                    group_meta=group_meta,
                )
                if rec:
                    records_buffer.append(rec)
                    total_records += 1
                    processed_groups += 1
                else:
                    skipped_groups += 1
                    self.logger.debug(
                        "Skipping PG group (%s, %s): no record built",
                        current_anchor,
                        current_run,
                    )
            except Exception as e:
                failed_groups += 1
                self.logger.debug(
                    "Failed PG group (%s, %s): %s",
                    current_anchor,
                    current_run,
                    e,
                )
            finally:
                current_rows = []

        with PgWriter(output_path, creator=creator, compression=self._compression) as writer:
            for batch in self._query_batched(grouped_sql, batch_size=query_batch_size):
                batch_df = batch.to_pandas()
                for row in batch_df.to_dict("records"):
                    # -- Validate group keys; skip null-like values ------
                    anchor_raw = row.get("anchor_protein")
                    if anchor_raw is None or (isinstance(anchor_raw, float) and pd.isna(anchor_raw)):
                        null_key_count += 1
                        continue
                    anchor = str(anchor_raw).strip()
                    if not anchor or anchor.lower() in ("none", "nan", "null"):
                        null_key_count += 1
                        continue

                    run_raw = row.get("run_file_name")
                    if run_raw is None or (isinstance(run_raw, float) and pd.isna(run_raw)):
                        null_key_count += 1
                        continue
                    run_file = str(run_raw).strip()
                    if not run_file or run_file.lower() in ("none", "nan", "null"):
                        null_key_count += 1
                        continue
                    # ----------------------------------------------------

                    if current_anchor is None:
                        current_anchor, current_run = anchor, run_file

                    if anchor != current_anchor or run_file != current_run:
                        _flush_current_group()
                        current_anchor, current_run = anchor, run_file

                    current_rows.append(row)

                    if len(records_buffer) >= write_batch_size:
                        self._track_scores(records_buffer)
                        writer.write_batch(records_buffer)
                        records_buffer.clear()

            _flush_current_group()
            if records_buffer:
                self._track_scores(records_buffer)
                writer.write_batch(records_buffer)

        bad_groups = skipped_groups + failed_groups
        summary = {
            "pg_groups_total": total_groups,
            "pg_groups_processed": processed_groups,
            "pg_groups_skipped": skipped_groups,
            "pg_groups_failed": failed_groups,
            "pg_bad_group_ratio": round((skipped_groups + failed_groups) / max(total_groups, 1), 4),
            "pg_null_anchor_keys": null_key_count,
        }
        self.logger.info("PG conversion summary: %s", summary)
        if bad_groups > 0:
            self.logger.warning(
                "PG conversion encountered problematic groups: %d/%d (skipped=%d, failed=%d)",
                bad_groups,
                total_groups,
                skipped_groups,
                failed_groups,
            )
        if null_key_count > 0:
            self.logger.warning(
                "PG conversion skipped %d feature rows with null-like anchor_protein or run_file_name keys",
                null_key_count,
            )
        if total_groups > 0 and processed_groups == 0:
            raise ValueError(
                "PG conversion failed: all protein groups were skipped or failed "
                f"(total={total_groups}, skipped={skipped_groups}, failed={failed_groups})."
            )
        if (
            total_groups >= min_groups_for_ratio_failure
            and total_groups > 0
            and (bad_groups / total_groups) > max_bad_group_ratio
        ):
            raise ValueError(
                "PG conversion failed: too many problematic protein groups "
                f"({bad_groups}/{total_groups}, skipped={skipped_groups}, failed={failed_groups}, "
                f"threshold={max_bad_group_ratio:.0%})."
            )

        self.logger.info(
            "PG conversion complete -> %s (%d records, streamed)",
            output_path,
            total_records,
        )

    # ------------------------------------------------------------------
    # Protein group metadata
    # ------------------------------------------------------------------

    def _build_protein_meta(
        self,
    ) -> tuple[dict[str, dict], dict[str, dict]]:
        """Build protein metadata from mzTab proteins table.

        Stores ``single_protein`` and ``indistinguishable_protein_group``
        rows keyed by their individual accession in *single_meta*.
        The ``indistinguishable_protein_group`` entries are **also** stored
        in *group_meta* by sorted, semicolon-joined accessions so that
        callers can look up the group entry directly.

        ``protein_details`` rows are skipped to match the dev version
        behaviour (they carry per-protein metadata but the dev SQL join
        never sees them).

        Returns:
            (single_meta, group_meta) tuple of dicts.
        """
        single_meta: dict[str, dict] = {}
        group_meta: dict[str, dict] = {}
        try:
            df = self._conn.execute("SELECT * FROM proteins").df()
            for row in df.to_dict("records"):
                result_type = str(row.get("opt_global_result_type", "single_protein"))
                accession = str(row.get("accession", ""))
                if not accession or accession == "null":
                    continue
                if result_type not in (
                    "single_protein",
                    "indistinguishable_protein_group",
                ):
                    continue

                entry = self._parse_protein_row(row, accession)

                if result_type == "single_protein":
                    if accession not in single_meta:
                        single_meta[accession] = entry
                else:
                    # indistinguishable_protein_group: store by individual
                    # accession (mzTab stores one accession per row, even
                    # for groups) so the lookup matches the dev SQL join.
                    if accession not in single_meta:
                        single_meta[accession] = entry
                    # Also keep the sorted joined key for group-level lookup
                    sorted_key = ";".join(sorted(self._normalize_pg_accessions(accession)))
                    if sorted_key:
                        group_meta[sorted_key] = entry
        except Exception as e:
            self.logger.error(f"Error loading protein metadata: {e}")

        return single_meta, group_meta

    def _parse_protein_row(self, row: dict, accession: str) -> dict:
        """Extract metadata fields from a single mzTab protein row."""
        description = row.get("description", "")
        description = description if pd.notna(description) and description != "null" else None

        # Gene names from description
        gg_accessions = []
        if description:
            gn = re.search(r"GN=([^\s]+)", str(description))
            if gn:
                gg_accessions = [gn.group(1)]

        # Global q-value
        global_qvalue = safe_float(
            row.get(
                "best_search_engine_score[1]",
                row.get("best_search_engine_score_1"),
            )
        )

        # Decoy
        is_decoy_raw = row.get(
            "opt_global_cv_pride:0000303_decoy_hit",
            row.get("opt_global_cv_PRIDE:0000303_decoy_hit", 0),
        )
        is_decoy = str(is_decoy_raw).strip() in ("1", "true", "True")

        # Sequence coverage
        seq_cov = safe_float(row.get("protein_coverage")) or 0.0

        # Sequence length (for iBAQ)
        seq_len = safe_float(row.get("sequence_length"))

        # Protein names
        pg_names = self._extract_protein_names(accession)

        return {
            "pg_names": pg_names,
            "gg_accessions": gg_accessions or None,
            "global_qvalue": global_qvalue,
            "is_decoy": is_decoy,
            "sequence_coverage": seq_cov,
            "sequence_length": seq_len,
        }

    # ------------------------------------------------------------------
    # Record building
    # ------------------------------------------------------------------

    def _build_single_pg(
        self,
        anchor_protein: str,
        run_file_name: str,
        features: pd.DataFrame,
        single_meta: dict[str, dict],
        group_meta: dict[str, dict],
    ) -> Optional[dict]:
        """Build a single PG record from aggregated features.

        Merges metadata from both single_protein and group rows:
        - q-value from single_protein (per-protein score)
        - decoy/coverage from whichever has them (prefer single_protein)
        """

        # Protein group accessions from first feature's pg_accessions.
        # Support list<dict>, list[str], semicolon-delimited strings, and null values.
        first_pg = features["pg_accessions"].iloc[0]
        pg_accessions = self._normalize_pg_accessions(first_pg)
        if not pg_accessions:
            self.logger.warning(
                "Empty/invalid pg_accessions for anchor=%s run=%s; falling back to anchor_protein",
                anchor_protein,
                run_file_name,
            )
            pg_accessions = [anchor_protein]

        # Lookup: anchor first, then each pg_accessions member.
        # The dev version's SQL join matches msstats ProteinName (which may
        # be semicolon-joined) to individual protein_groups entries.  We
        # replicate that by trying each member of pg_accessions.
        smeta = single_meta.get(anchor_protein, {})
        if not smeta:
            for acc in pg_accessions:
                smeta = single_meta.get(acc, {})
                if smeta:
                    break

        sorted_key = ";".join(sorted(pg_accessions))
        gmeta = group_meta.get(sorted_key, {})

        # Merge: prefer smeta (anchor or matched member), fall back to group
        pg_names_str = smeta.get("pg_names") or gmeta.get("pg_names")
        pg_names = pg_names_str.split(";") if pg_names_str else None

        gg_accessions = smeta.get("gg_accessions") or gmeta.get("gg_accessions")
        global_qvalue = smeta.get("global_qvalue") if smeta.get("global_qvalue") is not None else gmeta.get("global_qvalue")
        is_decoy = smeta.get("is_decoy", gmeta.get("is_decoy", False))
        seq_coverage = (
            smeta.get("sequence_coverage") if smeta.get("sequence_coverage") is not None else gmeta.get("sequence_coverage")
        )

        # Fallback: if neither dict had data, use feature-level info
        if not smeta and not gmeta:
            is_decoy_raw = features["is_decoy"].iloc[0]
            is_decoy = bool(is_decoy_raw) if pd.notna(is_decoy_raw) else False
            if "gg_accessions" in features.columns:
                first_gg = features["gg_accessions"].iloc[0]
                if isinstance(first_gg, list) and len(first_gg) > 0:
                    gg_accessions = first_gg

        # Aggregate intensities: sum per label across all features
        intensities = self._aggregate_intensities(features)

        # Peptide and feature counts
        unique_sequences = features["sequence"].nunique()
        total_sequences = len(features)

        feature_keys = set(zip(features["peptidoform"], features["charge"]))
        unique_features = len(feature_keys)
        total_features = len(features)

        # Peptides per protein
        peptides = [{"protein_name": acc, "peptide_count": unique_sequences} for acc in pg_accessions]

        # Additional scores
        additional_scores = []
        if seq_coverage is not None:
            additional_scores.append(
                {
                    "score_name": "sequence_coverage_percent",
                    "score_value": seq_coverage,
                    "higher_better": True,
                }
            )
        additional_scores.append(
            {
                "score_name": "peptide_count",
                "score_value": float(unique_sequences),
                "higher_better": True,
            }
        )

        return {
            "pg_accessions": pg_accessions,
            "pg_names": pg_names,
            "gg_accessions": (gg_accessions if isinstance(gg_accessions, list) else None),
            "gg_names": (
                gg_accessions if isinstance(gg_accessions, list) else None
            ),  # Gene symbols serve as both accession and name
            "gg_qvalue": None,
            "anchor_protein": anchor_protein,
            "run_file_name": run_file_name,
            "global_qvalue": global_qvalue,
            "pg_qvalue": global_qvalue,
            "intensities": intensities or None,
            "additional_intensities": self._compute_additional_intensities(
                features,
                intensities,
                smeta,
                gmeta,
                topn=self._topn,
                aggregation=self._aggregation,
            ),
            "is_decoy": is_decoy,
            "contaminant": False,
            "peptides": peptides,
            "peptide_counts": {
                "unique_sequences": unique_sequences,
                "total_sequences": total_sequences,
            },
            "feature_counts": {
                "unique_features": unique_features,
                "total_features": total_features,
            },
            "sequence_coverage": seq_coverage,
            "molecular_weight": None,
            "additional_scores": additional_scores or None,
            "cv_params": None,
        }

    def _aggregate_intensities(self, features: pd.DataFrame) -> list[dict]:
        """Sum feature intensities per label.

        Feature intensities are stored as list<intensity> where
        intensity = {label: string, intensity: float32}.
        """
        label_sums: dict[str, float] = {}
        for intensities in features["intensities"]:
            if intensities is None:
                continue
            for entry in intensities:
                label = entry.get("label", "LFQ")
                value = float(entry.get("intensity", 0.0) or 0.0)
                label_sums[label] = label_sums.get(label, 0.0) + value

        return [{"label": label, "intensity": total} for label, total in sorted(label_sums.items())]

    def _compute_additional_intensities(
        self,
        features: pd.DataFrame,
        intensities: list[dict],
        smeta: dict,
        gmeta: dict,
        topn: Optional[int] = None,
        aggregation: str = "sum",
    ) -> Optional[list[dict]]:
        """Compute TopN and iBAQ additional intensities per label.

        Returns schema-compatible ``list<additional_intensity>`` where
        each entry has ``{label, intensities: [{intensity_name, intensity_value}]}``.

        Args:
            topn: Number of top peptides.  ``None`` uses all peptides.
            aggregation: ``"sum"`` (default) or ``"mean"``.
        """
        # Flatten feature intensities to (peptidoform, label, intensity).
        # Use column arrays directly instead of iterrows() for speed.
        pf_col = features["peptidoform"].values
        ints_col = features["intensities"].values

        rows: list[tuple[str, str, float]] = []
        for idx in range(len(pf_col)):
            pf = str(pf_col[idx]) if pf_col[idx] is not None else ""
            ints = ints_col[idx]
            if ints is None:
                continue
            for entry in ints:
                label = entry.get("label", "LFQ")
                val = float(entry.get("intensity", 0.0) or 0.0)
                if val > 0:
                    rows.append((pf, label, val))

        if not rows:
            return None

        flat = pd.DataFrame(rows, columns=["peptidoform", "label", "intensity"])

        # Sequence length for iBAQ (prefer single_meta, fallback group_meta)
        seq_len = smeta.get("sequence_length") if smeta.get("sequence_length") else gmeta.get("sequence_length")

        # Determine intensity name
        intensity_name = f"Top{topn}" if topn is not None else "AllPeptides"

        result: list[dict] = []
        for label, grp in flat.groupby("label"):
            extra: list[dict] = []

            # Inline TopN: sort desc, optionally take top N, then aggregate
            sorted_ints = grp["intensity"].sort_values(ascending=False)
            if topn is not None:
                sorted_ints = sorted_ints.head(topn)

            if aggregation == "mean":
                topn_val = float(sorted_ints.mean())
            else:
                topn_val = float(sorted_ints.sum())

            extra.append({"intensity_name": intensity_name, "intensity_value": topn_val})

            # iBAQ = sum_intensity / sequence_length
            if seq_len and seq_len > 0:
                sum_int = float(grp["intensity"].sum())
                ibaq_val = sum_int / seq_len
                extra.append({"intensity_name": "ibaq", "intensity_value": ibaq_val})

            if extra:
                result.append({"label": str(label), "intensities": extra})

        return result or None

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _normalize_pg_accessions(raw_accessions) -> list[str]:
        """Normalize mixed pg_accessions shapes to a clean list of strings."""
        normalized: list[str] = []

        def _add(value) -> None:
            if value is None:
                return
            # Handle containers recursively
            if isinstance(value, (list, tuple)):
                for item in value:
                    _add(item)
                return
            if isinstance(value, dict):
                acc = value.get("accession") or value.get("pg_accession")
                if acc:
                    _add(acc)
                return
            # Scalar -- safe to use pd.isna
            try:
                if pd.isna(value):
                    return
            except (ValueError, TypeError):
                pass
            text = str(value).strip()
            if not text or text.lower() == "null":
                return
            for part in text.split(";"):
                part = part.strip()
                if part and part.lower() != "null":
                    normalized.append(part)

        items = list(raw_accessions) if isinstance(raw_accessions, (list, np.ndarray)) else [raw_accessions]
        for item in items:
            _add(item)
        return normalized

    @staticmethod
    def _extract_protein_names(accession: str) -> str:
        """Extract protein names from UniProt-style accession."""
        parts_list = []
        for acc in accession.split(";"):
            if "|" in acc:
                pieces = acc.split("|")
                if len(pieces) >= 3:
                    parts_list.append(pieces[2])
                elif len(pieces) >= 2:
                    parts_list.append(pieces[1])
                else:
                    parts_list.append(acc)
            else:
                parts_list.append(acc)
        return ";".join(parts_list)
