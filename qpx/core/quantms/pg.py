"""QuantMS-IO protein groups processing module."""

from __future__ import annotations

import logging
import os
import re
import time

import pandas as pd
import pyarrow as pa
from qpx.core.format import PG_SCHEMA
from qpx.core.quantms.mztab import MzTabIndexer
from qpx.utils.constants import MZTAB_PROTEIN_BEST_SEARCH_ENGINE_SCORE

from qpx.config import get_default_filters


class MzTabProteinGroups:
    """Handle protein groups in mzTab format with optimized quantification using composition pattern."""

    def __init__(self, mztab_indexer):
        """Initialize MzTabProteinGroups with an MzTabIndexer instance.

        Args:
            mztab_path: Path to mzTab file
        """
        self._indexer: MzTabIndexer = mztab_indexer

        # Initialize tracking lists first
        self._temp_files = []  # Track any temporary files created
        self._file_handles = []  # Track any open file handles
        # protein_details_df = self._indexer.get_protein_details()
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit with cleanup."""
        self.cleanup()

    def cleanup(self):
        """Clean up any temporary files and resources."""

        # Close any open file handles
        for file_handle in self._file_handles:
            try:
                if hasattr(file_handle, "close") and not file_handle.closed:
                    file_handle.close()
                    self.logger.debug("Closed file handle")
            except Exception as e:
                self.logger.warning(f"Error closing file handle: {e}")
        self._file_handles.clear()

        # Clean up any temporary files
        for temp_file in self._temp_files:
            try:
                if os.path.exists(temp_file):
                    os.unlink(temp_file)
                    self.logger.debug(f"Deleted temporary file: {temp_file}")
            except Exception as e:
                self.logger.warning(f"Error deleting temporary file {temp_file}: {e}")
        self._temp_files.clear()

    def __del__(self):
        """Destructor to ensure cleanup."""
        try:
            self.cleanup()
        except Exception as e:
            self.logger.warning(f"Exception during __del__ cleanup: {e}")

    def _is_protein_line(self, line_num: int) -> bool:
        """Check if a line is a protein data line."""
        try:
            # Use safe file opening with automatic cleanup
            with self._safe_file_open(self._mztab_path, "r") as file:
                for i, line in enumerate(file):
                    if i == line_num:
                        return line.startswith("PRT")
            return False
        except Exception:
            return False

    def _extract_protein_name(self, accession: str) -> str:
        """Extract protein name from accession."""
        if not accession or pd.isna(accession):
            return ""

        accessions = accession.split(";")
        acc_list = list()
        for acc in accessions:
            # Handle different accession formats
            if "|" in acc:
                # UniProt format: sp|P12345|PROT_HUMAN
                parts = acc.split("|")
                if len(parts) >= 3:
                    acc_list.append(parts[2])  # PROT_HUMAN
                elif len(parts) >= 2:
                    acc_list.append(parts[1])  # P12345
        if acc_list:
            return ";".join(acc_list)

        return accession

    def _extract_gene_names(self, description: str) -> list:
        """Extract gene names from protein description."""

        if not description or pd.isna(description):
            return []

        gene_names = []

        # Look for GN= pattern (common in UniProt)
        gn_match = re.search(r"GN=([^\s]+)", description)
        if gn_match:
            gene_names.append(gn_match.group(1))

        # Look for Gene_Symbol= pattern
        gs_match = re.search(r"Gene_Symbol=([^\s]+)", description)
        if gs_match:
            gene_names.append(gs_match.group(1))

        return list(set(gene_names))  # Remove duplicates

    def _safe_int_conversion(self, value, default=0):
        """Safely convert value to int with fallback."""
        try:
            if pd.isna(value) or value in ["null", "NULL", "", None]:
                return default
            return int(float(value))
        except (ValueError, TypeError):
            return default

    def _safe_file_open(self, file_path, mode="r"):
        """Safely open files with proper resource tracking."""
        import gzip
        from contextlib import contextmanager

        @contextmanager
        def managed_file():
            file_handle = None
            try:
                # Handle both compressed and uncompressed files
                if str(file_path).endswith(".gz"):
                    # For gzipped files, use gzip.open with text mode
                    file_handle = gzip.open(file_path, "rt" if "t" not in mode else mode, encoding="utf-8")
                else:
                    # For regular files
                    file_handle = open(file_path, mode, encoding="utf-8")

                # Track the file handle for cleanup (only if context manager doesn't close it)
                self._file_handles.append(file_handle)

                yield file_handle

            except Exception as e:
                self.logger.error(f"Error opening file {file_path}: {e}")
                raise
            finally:
                if file_handle:
                    # Remove from tracking list since we're closing it here
                    try:
                        self._file_handles.remove(file_handle)
                    except ValueError as e:
                        self.logger.warning(f"File handle already removed: {e}")

                    # Close the file
                    try:
                        if not file_handle.closed:
                            file_handle.close()
                    except Exception as e:
                        self.logger.warning(f"Exception closing file handle: {e}")

        return managed_file()

    def _convert_to_parquet_format(self, df: pd.DataFrame) -> pa.Table:
        """Convert DataFrame to parquet format using PG_SCHEMA."""
        if df.empty:
            # For empty DataFrames, create an empty table with the correct schema
            return pa.Table.from_arrays([pa.array([], type=field.type) for field in PG_SCHEMA], schema=PG_SCHEMA)
        else:
            return pa.Table.from_pandas(df, schema=PG_SCHEMA, preserve_index=False)

    def quantify_from_msstats_optimized(
        self,
        apply_filters: bool = True,
        filter_config: dict = None,
    ) -> pd.DataFrame:
        """Optimized protein quantification using DuckDB SQL aggregation.

        Args:
            apply_filters: Whether to apply QC filters from config.
            filter_config: Custom filter config dict. If None, uses default filters.
        """
        self.logger.info("[OPTIMIZED] Starting protein group quantification using DuckDB SQL")

        total_start = time.time()

        # Load filter configuration
        if apply_filters:
            filters = filter_config if filter_config else get_default_filters()
            self.logger.info(f"[FILTER] Using filter config: {filters.get('name', 'custom')}")
        else:
            filters = None

        # Step 1: Loading protein groups table from MzTabIndexer
        self.logger.info("[SETUP] Loading protein groups table from MzTabIndexer...")
        pg_start = time.time()
        protein_groups_info = self._load_protein_groups_table_optimized()
        pg_time = time.time() - pg_start
        self.logger.info(f"[SETUP] Created protein groups table with {len(protein_groups_info)} entries in {pg_time:.2f}s")

        # Step 2: Loading msstats data from MzTabIndexer for enhanced analysis
        self.logger.info("[DATA] Loading msstats data from MzTabIndexer...")
        msstats_start = time.time()

        # Step 3: Create joined view between msstats and protein groups
        self._create_msstats_protein_join_optimized(protein_groups_info)
        msstats_time = time.time() - msstats_start
        self.logger.info(f"[DATA] MzTabIndexer setup completed in {msstats_time:.2f}s")

        # Step 4: Process in batches using SQL aggregation
        self.logger.info("[PROCESS] Processing protein quantification in batches...")
        process_start = time.time()

        expanded_rows = []
        processed_files = 0

        for batch in self.get_sql_batchs():
            batch_start = time.time()

            try:
                batch_data = self.get_sql_batch_data(batch)

                self.logger.info(f"[SQL] Returned {len(batch_data)} rows for batch {batch}")

                if len(batch_data) > 0:
                    # Apply filters to batch data
                    if filters:
                        batch_data = self._apply_filters(batch_data, filters, self.logger)

                    if len(batch_data) > 0:
                        protein_row = self._create_optimized_protein_row(batch_data)
                        expanded_rows.extend(protein_row)

                    self.logger.info(f"[SUCCESS] Converted {len(batch_data)} SQL rows to protein rows")
                else:
                    self.logger.warning("[WARNING] No data returned from SQL")
                processed_files += len(batch_data)
                batch_time = time.time() - batch_start
                self.logger.info(f"[BATCH] Processed {len(batch_data)} rows in {batch_time:.2f}s ({processed_files} total)")

            except Exception as e:
                self.logger.error(f"[ERROR] SQL batch failed: {e}")
                self.logger.exception("Full traceback:")
                continue

        process_time = time.time() - process_start
        self.logger.info(f"[PROCESS] Completed quantification processing in {process_time:.2f}s")

        # Step 5: Convert to DataFrame
        self.logger.info("[CONVERT] Converting results to DataFrame...")
        df_start = time.time()

        if expanded_rows:
            result_df = pd.DataFrame(expanded_rows)
            self.logger.info(f"[CONVERT] Created DataFrame with {len(result_df)} rows and {len(result_df.columns)} columns")
        else:
            self.logger.warning("[WARNING] No data to convert - creating empty DataFrame")
            result_df = pd.DataFrame(
                columns=[
                    "pg_accessions",
                    "pg_names",
                    "gg_accessions",
                    "gg_names",
                    "reference_file_name",
                    "global_qvalue",
                    "pg_qvalue",
                    "intensities",
                    "additional_intensities",
                    "is_decoy",
                    "contaminant",
                    "peptides",
                    "anchor_protein",
                    "sequence_coverage",
                    "molecular_weight",
                    "additional_scores",
                    "cv_params",
                    "peptide_counts",
                    "feature_counts",
                ]
            )

        df_time = time.time() - df_start
        total_time = time.time() - total_start

        self.logger.info(f"[CONVERT] DataFrame conversion completed in {df_time:.2f}s")
        self.logger.info(f"[SUCCESS] OPTIMIZED quantification completed in {total_time:.2f}s total")

        # Clean up the temporary MzTabIndexer
        self._indexer.cleanup_duckdb()

        return result_df

    def _load_protein_groups_table_optimized(self):
        """Load protein groups lookup from mzTab protein section."""
        logger = self.logger

        protein_data = list()

        try:
            protein_df = self._indexer.get_proteins()
            total_rows = 0
            skipped_accession = 0
            skipped_type = 0
            result_types_seen = set()

            logger.info(f"[DEBUG] Protein DataFrame columns: {list(protein_df.columns)}")
            logger.info(f"[DEBUG] Protein DataFrame shape: {protein_df.shape}")

            for chunk in self.iter_in_chunks(protein_df):
                for _, row in chunk.iterrows():
                    total_rows += 1

                    result_type = row.get("opt_global_result_type", "single_protein")
                    result_types_seen.add(str(result_type))
                    accession = row.get("accession")

                    if pd.isna(accession) or not accession or accession == "null":
                        skipped_accession += 1
                        continue

                    # Skip other types
                    if result_type not in [
                        "single_protein",
                        "indistinguishable_protein_group",
                    ]:
                        skipped_type += 1
                        continue

                    pg_accessions = accession.strip()
                    anchor_protein = pg_accessions
                    description = self.check_clean_null(row["description"], None)

                    protein_data.append(
                        {
                            "result_type": result_type,
                            "pg_accessions": pg_accessions,
                            "anchor_protein": anchor_protein,
                            "description": description,
                            "sequence_coverage": self.check_clean_null(row["protein_coverage"], 0),
                            "global_qvalue": self.check_clean_null(row[MZTAB_PROTEIN_BEST_SEARCH_ENGINE_SCORE], 0),
                            "peptide_count": self.check_clean_null(row["opt_global_nr_found_peptides"], 0),
                            "is_decoy": self._safe_int_conversion(
                                self.check_clean_null(row.get("opt_global_cv_PRIDE:0000303_decoy_hit"), 0)
                            ),
                            "sequence_length": self.check_clean_null(row.get("sequence_length"), 0),
                            "pg_names": self._extract_protein_name(pg_accessions),
                            "gg_accessions": self._extract_gene_names(description),
                        }
                    )

        except Exception as e:
            logger.error(f"Error loading protein groups table: {e}")
            logger.exception("Full traceback:")

        logger.info(f"[DEBUG] Result types seen: {result_types_seen}")
        logger.info(
            f"[DEBUG] Skipped {skipped_accession} rows due to missing accession, {skipped_type} rows due to result_type filter"
        )
        logger.info(f"Loaded {len(protein_data)} protein groups from {total_rows} total rows")
        return protein_data

    def check_clean_null(self, value, fill_value):
        return value if pd.notna(value) and value != "null" else fill_value

    def iter_in_chunks(self, df, chunksize=100000):
        for start in range(0, len(df), chunksize):
            yield df.iloc[start : start + chunksize]

    def _create_msstats_protein_join_optimized(self, protein_groups_info):
        """Create optimized join view between msstats and protein groups in DuckDB."""
        logger = self.logger

        try:
            if protein_groups_info:
                protein_df = pd.DataFrame(protein_groups_info)

                self._indexer._duckdb.execute("DROP TABLE IF EXISTS protein_groups")
                self._indexer._duckdb.from_df(protein_df).create("protein_groups")

                # Create index for better join performance
                self._indexer._duckdb.execute("CREATE INDEX IF NOT EXISTS idx_protein_name ON protein_groups(pg_accessions)")

                # Get 'unique' from PSM table
                unique_peptide_df = self._indexer.get_unique_from_psm_table()

                self._indexer._duckdb.execute("DROP TABLE IF EXISTS unique_peptide")
                self._indexer._duckdb.from_df(unique_peptide_df).create("unique_peptide")

                # Simplified join - use exact match first, then fallback for unmatched
                self._indexer._duckdb.execute("""
                    DROP VIEW IF EXISTS processed_msstats_with_pg;
                    CREATE VIEW processed_msstats_with_pg AS
                    SELECT
                        m.*,
                        COALESCE(pg.anchor_protein, m.pg_accessions) AS anchor_protein,
                        pg.result_type,
                        pg.description,
                        pg.sequence_coverage,
                        pg.global_qvalue,
                        pg.peptide_count,
                        pg.is_decoy,
                        pg.sequence_length,
                        pg.pg_names,
                        pg.gg_accessions,
                        unpep."unique"
                    FROM msstats m
                    LEFT JOIN protein_groups pg ON m.pg_accessions = pg.pg_accessions
                    LEFT JOIN unique_peptide unpep ON m.pg_accessions = unpep.pg_accessions AND m.peptidoform = unpep.peptidoform
                """)

                # Log statistics
                count = self._indexer._duckdb.execute("SELECT COUNT(*) FROM msstats").fetchone()[0]
                matched_count = self._indexer._duckdb.execute(
                    "SELECT COUNT(*) FROM processed_msstats_with_pg WHERE anchor_protein IS NOT NULL"
                ).fetchone()[0]
                logger.info(f"[DATA] Created msstats view with {count} rows, {matched_count} matched to protein groups")

            else:
                logger.warning("[WARNING] No protein groups data to create lookup table")

        except Exception as e:
            logger.error(f"Error creating msstats protein join: {e}")
            raise

    def _apply_filters(self, batch_data: pd.DataFrame, filters: dict, logger) -> pd.DataFrame:
        """Apply QC filters to batch data based on filter configuration.

        Args:
            batch_data: DataFrame with peptide/protein data
            filters: Filter configuration dict
            logger: Logger instance

        Returns:
            Filtered DataFrame
        """
        original_count = len(batch_data)

        # Intensity filters
        intensity_filters = filters.get("intensity", {})
        if intensity_filters.get("remove_zero_intensity", False):
            batch_data = batch_data[batch_data["intensity"] > 0]

        # Peptide filters
        peptide_filters = filters.get("peptide", {})
        min_len = peptide_filters.get("min_peptide_length")
        max_len = peptide_filters.get("max_peptide_length")
        if min_len and "peptidoform" in batch_data.columns:
            batch_data = batch_data[batch_data["peptidoform"].str.len() >= min_len]
        if max_len and "peptidoform" in batch_data.columns:
            batch_data = batch_data[batch_data["peptidoform"].str.len() <= max_len]

        # Protein filters
        protein_filters = filters.get("protein", {})
        contaminant_patterns = protein_filters.get("contaminant_patterns", [])
        if protein_filters.get("remove_contaminants", False) and contaminant_patterns:
            pattern = "|".join(contaminant_patterns)
            batch_data = batch_data[~batch_data["anchor_protein"].str.contains(pattern, case=False, na=False)]

        filtered_count = len(batch_data)
        if filters.get("log_filtered_counts", False) and original_count != filtered_count:
            logger.info(f"[FILTER] Filtered {original_count - filtered_count} rows ({original_count} -> {filtered_count})")

        return batch_data

    @staticmethod
    def _prepare_peptide_data(channel_group, sample_accession):
        """Prepare peptide data for quantification."""
        peptide_data = channel_group[["pg_accessions", "peptidoform", "intensity"]].copy()
        peptide_data = peptide_data.rename(
            columns={
                "pg_accessions": "ProteinName",
                "peptidoform": "PeptideSequence",
                "intensity": "NormIntensity",
            }
        )
        peptide_data["SampleID"] = sample_accession
        return peptide_data

    @staticmethod
    def _compute_additional_scores(group, max_intensity, avg_intensity, unique_peptide_count):
        """Compute additional quality scores for a protein group."""
        extra_scores = []
        sequence_coverage = float(group["sequence_coverage"].iloc[0])
        if sequence_coverage != 0:
            extra_scores.append(
                {
                    "score_name": "sequence_coverage_percent",
                    "score_value": sequence_coverage,
                    "higher_better": True,
                }
            )
        extra_scores.append(
            {
                "score_name": "peptide_count",
                "score_value": float(
                    max(
                        float(group["peptide_count"].iloc[0]),
                        unique_peptide_count,
                    )
                ),
                "higher_better": True,
            }
        )
        extra_scores.append(
            {
                "score_name": "max_intensity",
                "score_value": max_intensity,
                "higher_better": True,
            }
        )
        extra_scores.append(
            {
                "score_name": "avg_intensity",
                "score_value": avg_intensity,
                "higher_better": True,
            }
        )
        return extra_scores

    def _create_optimized_protein_row(self, batch_data):
        result = []
        for (anchor_protein, reference_file_name), group in batch_data.groupby(["anchor_protein", "reference_file_name"]):
            unique_group = group[group["unique"].isin([1, "1"])]
            peptide_count = {
                "unique_sequences": unique_group["peptidoform"].nunique(),
                "total_sequences": group["peptidoform"].nunique(),
            }
            feature_count = {
                "unique_features": len(set(zip(unique_group["peptidoform"], unique_group["charge"]))),
                "total_features": len(set(zip(group["peptidoform"], group["charge"]))),
            }

            if float(group["peptide_count"].iloc[0]) != 0:
                unique_peptide_count = float(group["peptide_count"].iloc[0])
            else:
                unique_peptide_count = group["peptidoform"].nunique()

            max_intensity = group["intensity"].max()
            avg_intensity = group["intensity"].mean()

            # Peptide counts per protein in the group (consistent with other converters)
            pg_acc_raw = group["pg_accessions"].iloc[0]
            pg_acc_list = pg_acc_raw.split(";") if isinstance(pg_acc_raw, str) else [anchor_protein]
            peptides = [{"protein_name": acc.strip(), "peptide_count": unique_peptide_count} for acc in pg_acc_list]

            intensities = []
            additional_intensities = []

            for channel, channel_group in group.groupby("channel"):
                # sample_accession
                sample_accessions = channel_group["sample_accession"].unique()
                if len(sample_accessions) != 1:
                    raise ValueError("number of sample_accession is more than 1")
                sample_accession = sample_accessions[0]

                # intensity sum
                sum_intensity = float(channel_group["intensity"].sum())

                intensities.append(
                    {
                        "sample_accession": sample_accession,
                        "channel": channel,
                        "intensity": sum_intensity,
                    }
                )

                extra_scores = self._compute_additional_scores(group, max_intensity, avg_intensity, unique_peptide_count)

            # Extract gg_accessions as list
            gg_acc_raw = group["gg_accessions"].iloc[0]
            gg_accessions_list = None
            if gg_acc_raw is not None:
                if isinstance(gg_acc_raw, str) and gg_acc_raw:
                    gg_accessions_list = gg_acc_raw.split(";")
                elif hasattr(gg_acc_raw, "__len__") and len(gg_acc_raw) > 0:
                    gg_accessions_list = list(gg_acc_raw)

            # Extract sequence_coverage (already in extra_scores, but also as separate field)
            seq_coverage = group["sequence_coverage"].iloc[0]
            seq_coverage_val = float(seq_coverage) if pd.notna(seq_coverage) else None

            result.append(
                {
                    "pg_accessions": group["pg_accessions"].iloc[0].split(";"),
                    "pg_names": group["pg_names"].iloc[0].split(";"),
                    "gg_accessions": gg_accessions_list,
                    "gg_names": None,  # Not available in mzTab
                    "reference_file_name": reference_file_name,
                    "global_qvalue": float(group["global_qvalue"].iloc[0]),
                    "pg_qvalue": None,  # Not available in mzTab (DIA-NN specific)
                    "intensities": intensities,
                    "additional_intensities": additional_intensities,
                    "is_decoy": int(group["is_decoy"].iloc[0]),
                    "contaminant": 0,  # mzTab doesn't have contaminant info
                    "peptides": peptides,
                    "anchor_protein": anchor_protein.split(";")[0],
                    "sequence_coverage": seq_coverage_val,
                    "molecular_weight": None,  # Not available in mzTab
                    "additional_scores": extra_scores,
                    "cv_params": None,  # Not available in mzTab
                    "peptide_counts": peptide_count,
                    "feature_counts": feature_count,
                }
            )

        return result

    def get_sql_batchs(self, file_num: int = 2):
        try:
            query = """
                SELECT DISTINCT reference_file_name
                FROM processed_msstats_with_pg
                ORDER BY reference_file_name
            """
            files_df = self._indexer._duckdb.execute(query).df()
            unique_files = files_df["reference_file_name"].tolist()

            for i in range(0, len(unique_files), file_num):
                yield unique_files[i : i + file_num]

        except Exception as e:
            self.logger.error(f"Error getting file batches: {e}")
            return

    def get_sql_batch_data(self, file_batch):
        batch_query = """
            SELECT
                anchor_protein,
                pg_accessions,
                result_type,
                description,
                sequence_coverage,
                global_qvalue,
                peptide_count,
                is_decoy,
                sequence_length,
                reference_file_name,
                channel,
                pg_names,
                gg_accessions,
                peptidoform,
                intensity,
                sample_accession,
                charge,
                "unique"
            FROM processed_msstats_with_pg
            WHERE reference_file_name = ANY($1)
            AND anchor_protein IS NOT NULL
            AND intensity > 0
        """
        batch_data = self._indexer._duckdb.execute(batch_query, [file_batch]).df()

        return batch_data
