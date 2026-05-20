"""QuantMS Feature adapter -- mzTab + MSstats to feature.parquet.

Loads an mzTab file and its companion MSstats input into DuckDB, joins
them via SQL, and streams the result through ``FeatureWriter``.

Key schema changes:
    - ``reference_file_name`` -> ``run_file_name``
    - ``precursor_charge`` (int32) -> ``charge`` (int16)
    - ``is_decoy`` (int32, 0/1) -> ``is_decoy`` (bool)
    - ``scan`` (string) -> ``scan`` (list<int32>)
    - ``start_ion_mobility``/``stop_ion_mobility`` -> ``ion_mobility_start``/``ion_mobility_stop``
    - Intensities struct: ``{sample_accession, channel, intensity}`` ->
      ``{label, intensity}``
"""

from __future__ import annotations

import json
import logging
import re
from typing import Optional

import pandas as pd

from qpx.converters.base import BaseConverter, resolve_columns
from qpx.converters.mappings import get_field_mappings
from qpx.converters.mztab import (
    extract_modifications,
    extract_ms_runs,
    load_msstats,
    load_mztab_sections,
)
from qpx.converters.ptm import from_proforma
from qpx.converters.utils import (
    parse_scan_numbers,
    resolve_run_file,
    safe_float,
)
from qpx.core.cleavage import count_missed_cleavages
from qpx.core.sql import sql_build, validate_identifier, validate_table
from qpx.writers.feature import FeatureWriter

logger = logging.getLogger(__name__)

# Derive field map from constants
_FEATURE_MAP = get_field_mappings("quantms", "feature")


class QuantmsFeatureAdapter(BaseConverter):
    """Convert QuantMS mzTab + MSstats data to ``feature.parquet``.

    Usage::

        with QuantmsFeatureAdapter() as adapter:
            adapter.convert(
                mztab_path="results.mzTab",
                msstats_path="msstats_in.csv",
                output_path="feature.parquet",
            )
    """

    def __init__(self, **kwargs):
        """Initialize adapter with an empty peptide-protein map."""
        super().__init__(**kwargs)
        self._pep_protein_map: dict[str, str] = {}

    def convert(
        self,
        mztab_path: str,
        msstats_path: str,
        output_path: str,
        sdrf_path: Optional[str] = None,
        file_batch_size: int = 10,
        chunksize: int = 500_000,
        creator: str = "quantms",
        enzyme_name: str | None = None,
    ) -> None:
        """Run the mzTab+MSstats -> feature.parquet conversion.

        Args:
            mztab_path: Path to the mzTab file.
            msstats_path: Path to the MSstats input file.
            output_path: Destination Parquet path.
            sdrf_path: Optional SDRF path for sample mapping.
            file_batch_size: Number of run files to batch.
            chunksize: Rows per streaming batch.
            creator: Creator tag in Parquet metadata.
            enzyme_name: Enzyme name for computing missed cleavages.
        """
        self._enzyme_name = enzyme_name
        self._pep_protein_map: dict[str, str] = {}
        # Step 1: Load data into DuckDB (skip if already loaded)
        if not self._table_exists("psms"):
            load_mztab_sections(self._conn, mztab_path)
        if not self._table_exists("msstats"):
            load_msstats(self._conn, msstats_path)

        # Step 2: Extract auxiliary lookups
        ms_runs = extract_ms_runs(self._conn)
        self._modifications_meta = extract_modifications(self._conn)

        # Step 3: Determine experiment type (LFQ or TMT)
        experiment_type = self._detect_experiment_type()

        # Step 4: For LFQ, use optimized SQL-first path
        if experiment_type == "LFQ":
            self._convert_lfq_fast(
                output_path,
                ms_runs,
                creator,
                chunksize,
            )
        else:
            # Isobaric path uses the original dict-based approach
            psm_lookup = self._build_psm_lookup(ms_runs)
            protein_qvalue_map = self._build_protein_qvalue_map()
            protein_gene_map = self._build_protein_gene_map()

            self.logger.info("Aggregating MSstats data and writing features ...")
            with FeatureWriter(output_path, creator=creator, compression=self._compression) as writer:
                for batch_df in self._iter_feature_batches(file_batch_size):
                    records = self._transform_batch(
                        batch_df,
                        psm_lookup,
                        protein_qvalue_map,
                        protein_gene_map,
                        experiment_type,
                        ms_runs,
                    )
                    if records:
                        self._track_scores(records)
                        writer.write_batch(records)

        self.logger.info(f"Feature conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Fast LFQ path: SQL-first with DuckDB JOINs
    # ------------------------------------------------------------------

    def _convert_lfq_fast(
        self,
        output_path: str,
        ms_runs: dict[int, str],
        creator: str,
        chunksize: int,
    ) -> None:
        """Optimized LFQ feature conversion using DuckDB SQL JOINs.

        Instead of pulling rows into Python and looping, this method:
        1. Loads PSM lookup, protein q-values, gene maps as DuckDB tables
        2. Builds a ProForma lookup table (unique peptidoforms only)
        3. JOINs everything in SQL — string ops, mass error, etc. computed in SQL
        4. Streams pre-computed rows; Python only assembles nested structs
        """
        self.logger.info("Using fast SQL-first LFQ path ...")

        # Resolve MSstats column names
        msstats_cols = {
            c[0]
            for c in self._conn.execute(
                "SELECT column_name FROM information_schema.columns WHERE table_name='msstats'"
            ).fetchall()
        }
        col_map = resolve_columns(_FEATURE_MAP, msstats_cols)
        pf_col = col_map.get("peptidoform", "PeptideSequence")
        prot_col = col_map.get("pg_accessions", "ProteinName")
        ref_col = col_map.get("run_file_name", "Reference")
        charge_col = col_map.get("charge", "Charge")
        intensity_col = col_map.get("intensity", "Intensity")
        rt_col = col_map.get("rt", "RetentionTime")
        has_rt = rt_col in msstats_cols

        # Step 1: Load PSM lookup as DuckDB table
        self._load_psm_lookup_table(ms_runs)

        # Step 2: Load protein q-value and gene maps as DuckDB tables
        self._load_protein_tables()

        # Step 3: Build ProForma lookup (distinct peptidoforms only)
        self._load_proforma_lookup(pf_col)

        # Step 3b: Build peptide → protein lookup from mzTab PEP section
        self._pep_protein_map = self._build_peptide_protein_map()

        # Step 4: Build the big JOIN query
        # NOTE: Column names (pf_col, ref_col, etc.) come from the internal
        # _FEATURE_MAP dict resolved against actual DuckDB column names — they
        # are NOT user-supplied, so f-string interpolation is safe here.
        rt_expr = (
            sql_build(
                "TRY_CAST(m.$rt_col AS DOUBLE)",
                rt_col=validate_identifier(rt_col),
            )
            if has_rt
            else "NULL"
        )

        q_pf = validate_identifier(pf_col)
        q_ref = validate_identifier(ref_col)
        q_chg = validate_identifier(charge_col)
        q_int = validate_identifier(intensity_col)
        q_prot = validate_identifier(prot_col)

        sql = sql_build(
            """
            SELECT
                pf.peptidoform AS peptidoform,
                regexp_replace(
                    regexp_replace(upper(CAST(pf.peptidoform AS VARCHAR)), '\\[.*?\\]', '', 'g'),
                    '[^A-Z]',
                    '',
                    'g'
                ) AS sequence,
                split_part(CAST(m.$ref_col AS VARCHAR), '.', 1) AS run_file_name,
                COALESCE(TRY_CAST(m.$chg_col AS INTEGER), 0) AS charge,
                COALESCE(TRY_CAST(m.$int_col AS DOUBLE), 0.0) AS intensity,
                $rt_expr AS msstats_rt,
                m.$prot_col AS protein_name,
                -- PSM lookup fields
                p.pep AS psm_pep,
                p.calc_mz AS psm_calc_mz,
                p.obs_mz AS psm_obs_mz,
                p.is_decoy AS psm_is_decoy,
                p.scan AS psm_scan,
                p.id_run_file_name AS psm_id_run,
                p.rt AS psm_rt,
                -- Protein q-value
                pq.qvalue AS pg_global_qvalue,
                -- Gene names
                pg.gene_name AS gene_name,
                -- ProForma modifications (JSON)
                pf.modifications_json AS modifications_json
            FROM msstats m
            LEFT JOIN _psm_lookup p
                ON split_part(CAST(m.$ref_col AS VARCHAR), '.', 1) = p.run_file_name
                AND CAST(m.$pf_col AS VARCHAR) = p.peptidoform
                AND CAST(COALESCE(TRY_CAST(m.$chg_col AS INTEGER), 0) AS VARCHAR) = p.charge
            LEFT JOIN _protein_qvalues pq
                ON split_part(CAST(m.$prot_col AS VARCHAR), ';', 1) = pq.accession
            LEFT JOIN _protein_genes pg
                ON split_part(CAST(m.$prot_col AS VARCHAR), ';', 1) = pg.accession
            LEFT JOIN _proforma_lookup pf
                ON CAST(m.$pf_col AS VARCHAR) = pf.raw_peptidoform
            """,
            pf_col=q_pf,
            ref_col=q_ref,
            chg_col=q_chg,
            int_col=q_int,
            prot_col=q_prot,
            rt_expr=rt_expr,
        )

        # Step 5: Stream results and build feature records
        self.logger.info("Streaming SQL-joined feature rows ...")
        total = self._conn.execute("SELECT COUNT(*) FROM msstats").fetchone()[0]
        self.logger.info("Total MSstats rows to process: %d", total)

        with FeatureWriter(output_path, creator=creator, compression=self._compression) as writer:
            offset = 0
            batch_num = 0
            while offset < total:
                batch_sql = sql_build(
                    "$base LIMIT $lim OFFSET $off",
                    base=sql,
                    lim=str(int(chunksize)),
                    off=str(int(offset)),
                )
                rows = self._conn.execute(batch_sql).fetchall()
                if not rows:
                    break

                records = self._rows_to_feature_records(rows)
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)

                batch_num += 1
                offset += chunksize
                if batch_num % 5 == 0:
                    self.logger.info(
                        "Processed %d / %d rows (%.1f%%)",
                        min(offset, total),
                        total,
                        100 * min(offset, total) / total,
                    )

        # Clean up temp tables
        for t in ("_psm_lookup", "_protein_qvalues", "_protein_genes", "_proforma_lookup"):
            self._conn.execute(sql_build("DROP TABLE IF EXISTS $t", t=validate_table(t)))

    def _load_psm_lookup_table(self, ms_runs: dict[int, str]) -> None:
        """Load PSM data as a DuckDB table for SQL JOINs."""
        actual_cols = {
            c[0]
            for c in self._conn.execute("SELECT column_name FROM information_schema.columns WHERE table_name='psms'").fetchall()
        }

        # Peptidoform column
        pf_lo = "opt_global_cv_ms:1000889_peptidoform_sequence"
        pf_hi = "opt_global_cv_MS:1000889_peptidoform_sequence"
        pf_col = pf_lo if pf_lo in actual_cols else (pf_hi if pf_hi in actual_cols else None)

        # Decoy column
        dec_lo = "opt_global_cv_ms:1002217_decoy_peptide"
        dec_hi = "opt_global_cv_MS:1002217_decoy_peptide"
        dec_col = dec_lo if dec_lo in actual_cols else (dec_hi if dec_hi in actual_cols else None)

        # PEP column
        pep_col = None
        for candidate in ["opt_global_posterior_error_probability_score", "opt_global_Posterior_Error_Probability_score"]:
            if candidate in actual_cols:
                pep_col = candidate
                break

        pf_expr = (
            sql_build(
                "CAST($col AS VARCHAR)",
                col=validate_identifier(pf_col),
            )
            if pf_col
            else "''"
        )
        dec_expr = (
            sql_build(
                "CAST($col AS VARCHAR)",
                col=validate_identifier(dec_col),
            )
            if dec_col
            else "'0'"
        )
        pep_expr = (
            sql_build(
                "TRY_CAST($col AS DOUBLE)",
                col=validate_identifier(pep_col),
            )
            if pep_col
            else "NULL"
        )
        rt_expr = "TRY_CAST(retention_time AS DOUBLE)" if "retention_time" in actual_cols else "NULL"

        # Build ms_run mapping table
        self._conn.execute("""
            CREATE OR REPLACE TABLE _ms_runs (idx INTEGER, file_stem VARCHAR)
        """)
        if ms_runs:
            for idx, stem in ms_runs.items():
                self._conn.execute(
                    "INSERT INTO _ms_runs VALUES (?, ?)",
                    [idx, stem],
                )

        # Build PSM lookup table with first occurrence per (run, peptidoform, charge)
        self._conn.execute(
            sql_build(
                """
            CREATE OR REPLACE TABLE _psm_lookup AS
            WITH psm_enriched AS (
                SELECT
                    spectra_ref,
                    $pf_expr AS peptidoform,
                    CAST(charge AS VARCHAR) AS charge,
                    $pep_expr AS pep,
                    TRY_CAST(calc_mass_to_charge AS DOUBLE) AS calc_mz,
                    TRY_CAST(exp_mass_to_charge AS DOUBLE) AS obs_mz,
                    $dec_expr AS is_decoy_raw,
                    CAST(accession AS VARCHAR) AS accession,
                    $rt_expr AS rt,
                    -- Extract ms_run index from spectra_ref
                    TRY_CAST(
                        regexp_extract(CAST(spectra_ref AS VARCHAR), '\\[(\\d+)\\]', 1)
                        AS INTEGER
                    ) AS ms_run_idx,
                    -- Extract scan number
                    TRY_CAST(
                        regexp_extract(CAST(spectra_ref AS VARCHAR), '(?:scan|index|spectrum)=(\\d+)', 1)
                        AS INTEGER
                    ) AS scan_number,
                    ROW_NUMBER() OVER (
                        PARTITION BY
                            r.file_stem,
                            $pf_expr,
                            CAST(charge AS VARCHAR)
                        ORDER BY spectra_ref
                    ) AS rn
                FROM psms
                LEFT JOIN _ms_runs r ON TRY_CAST(
                    regexp_extract(CAST(spectra_ref AS VARCHAR), '\\[(\\d+)\\]', 1)
                    AS INTEGER
                ) = r.idx
            )
            SELECT
                split_part(COALESCE(r.file_stem, ''), '.', 1) AS run_file_name,
                peptidoform,
                charge,
                pep,
                calc_mz,
                obs_mz,
                CASE WHEN TRIM(is_decoy_raw) = '1' THEN TRUE ELSE FALSE END AS is_decoy,
                accession,
                scan_number AS scan,
                split_part(COALESCE(r.file_stem, ''), '.', 1) AS id_run_file_name,
                rt
            FROM psm_enriched pe
            LEFT JOIN _ms_runs r ON pe.ms_run_idx = r.idx
            WHERE rn = 1
            """,
                pf_expr=pf_expr,
                dec_expr=dec_expr,
                pep_expr=pep_expr,
                rt_expr=rt_expr,
            )
        )

        count = self._conn.execute("SELECT COUNT(*) FROM _psm_lookup").fetchone()[0]
        self.logger.info("PSM lookup table: %d entries", count)

        # Clean up ms_runs temp table
        self._conn.execute("DROP TABLE IF EXISTS _ms_runs")

    def _load_protein_tables(self) -> None:
        """Load protein q-value and gene maps as DuckDB tables."""
        # Protein q-values
        cols = {
            c[0]
            for c in self._conn.execute(
                "SELECT column_name FROM information_schema.columns WHERE table_name='proteins'"
            ).fetchall()
        }
        score_col = None
        for candidate in ["best_search_engine_score[1]", "best_search_engine_score_1"]:
            if candidate in cols:
                score_col = candidate
                break

        if score_col:
            q_score = validate_identifier(score_col)
            self._conn.execute(
                sql_build(
                    """CREATE OR REPLACE TABLE _protein_qvalues AS
                SELECT
                    CAST(accession AS VARCHAR) AS accession,
                    TRY_CAST($sc AS DOUBLE) AS qvalue
                FROM proteins
                WHERE accession IS NOT NULL AND $sc IS NOT NULL""",
                    sc=q_score,
                )
            )
        else:
            self._conn.execute("CREATE OR REPLACE TABLE _protein_qvalues (accession VARCHAR, qvalue DOUBLE)")

        # Protein -> gene map
        if "description" in cols:
            self._conn.execute("""
                CREATE OR REPLACE TABLE _protein_genes AS
                SELECT
                    CAST(accession AS VARCHAR) AS accession,
                    regexp_extract(CAST(description AS VARCHAR), 'GN=([^\\s]+)', 1) AS gene_name
                FROM proteins
                WHERE accession IS NOT NULL
                  AND description IS NOT NULL
                  AND CAST(description AS VARCHAR) != 'null'
                  AND regexp_extract(CAST(description AS VARCHAR), 'GN=([^\\s]+)', 1) != ''
            """)
        else:
            self._conn.execute("CREATE OR REPLACE TABLE _protein_genes (accession VARCHAR, gene_name VARCHAR)")

        qv_count = self._conn.execute("SELECT COUNT(*) FROM _protein_qvalues").fetchone()[0]
        gene_count = self._conn.execute("SELECT COUNT(*) FROM _protein_genes").fetchone()[0]
        self.logger.info("Protein tables: %d q-values, %d gene mappings", qv_count, gene_count)

    def _load_proforma_lookup(self, pf_col: str) -> None:
        """Extract unique peptidoforms, parse ProForma once, load as DuckDB table.

        For 7.6M MSstats rows there may be only ~200K unique peptidoforms,
        so parsing 200K instead of 7.6M is a huge win.
        """
        # Get distinct peptidoforms from MSstats
        q_pf = validate_identifier(pf_col)
        distinct = self._conn.execute(
            sql_build(
                """SELECT DISTINCT CAST($pf AS VARCHAR) AS peptidoform
            FROM msstats
            WHERE $pf IS NOT NULL""",
                pf=q_pf,
            )
        ).fetchall()

        self.logger.info("Parsing %d unique peptidoforms ...", len(distinct))

        mods_meta = self._modifications_meta
        records = []
        for (peptidoform,) in distinct:
            if not peptidoform or peptidoform == "null":
                continue
            sequence = re.sub(r"[^A-Z]", "", peptidoform.upper())
            if peptidoform != sequence:
                peptidoform_profoma, mods = from_proforma(
                    peptidoform,
                    sequence,
                    meta=mods_meta,
                )
                mods_json = json.dumps(mods) if mods else None
            else:
                mods_json = None
                peptidoform_profoma = peptidoform
            records.append((peptidoform, peptidoform_profoma, mods_json))

        # Load into DuckDB
        if records:
            import pandas as _pd

            df = _pd.DataFrame(
                records,
                columns=["raw_peptidoform", "peptidoform", "modifications_json"],
            )
            self._conn.execute("DROP TABLE IF EXISTS _proforma_lookup")
            self._conn.from_df(df).create("_proforma_lookup")
        else:
            self._conn.execute("""
            CREATE OR REPLACE TABLE _proforma_lookup (
                raw_peptidoform VARCHAR,
                peptidoform VARCHAR,
                modifications_json VARCHAR
            )
            """)
        self.logger.info("ProForma lookup table: %d entries", len(records))

    def _build_peptide_protein_map(self) -> dict[str, str]:
        """Build peptide sequence → single protein accession map from mzTab PEP section.

        The mzTab PEP section contains the protein inference result: each
        peptide is assigned to exactly one protein accession (including razor
        peptide resolution).  This map is used in ``_rows_to_feature_records``
        to resolve protein groups (``A;B``) to the correct single accession.

        Only unambiguous peptides (single protein) are included.

        Returns
        -------
        dict[str, str]
            Dict mapping plain sequence (uppercase, letters only) to single
            protein accession.

        """
        if not self._table_exists("peptides"):
            self.logger.info("No mzTab peptides table — skipping peptide protein map")
            return {}
        pep_cols = {
            c[0]
            for c in self._conn.execute(
                "SELECT column_name FROM information_schema.columns WHERE table_name='peptides'"
            ).fetchall()
        }
        if "sequence" not in pep_cols or "accession" not in pep_cols:
            self.logger.info("Peptides table lacks sequence/accession — skipping")
            return {}
        rows = self._conn.execute("""
            WITH pep_resolved AS (
                SELECT
                    regexp_replace(upper(CAST(sequence AS VARCHAR)), '[^A-Z]', '', 'g')
                        AS pep_sequence,
                    CAST(accession AS VARCHAR) AS pep_accession
                FROM peptides
                WHERE accession IS NOT NULL
                  AND CAST(accession AS VARCHAR) != 'null'
                  AND sequence IS NOT NULL
            )
            SELECT pep_sequence, MIN(pep_accession) AS pep_accession
            FROM pep_resolved
            GROUP BY pep_sequence
            HAVING COUNT(DISTINCT pep_accession) = 1
        """).fetchall()

        pep_map = dict(rows)
        self.logger.info("Peptide-protein map: %d entries (unambiguous)", len(pep_map))
        return pep_map

    def _rows_to_feature_records(self, rows: list[tuple]) -> list[dict]:
        """Convert pre-joined SQL rows to feature record dicts.

        Column order from the SQL query:
        0: peptidoform, 1: sequence, 2: run_file_name, 3: charge,
        4: intensity, 5: msstats_rt, 6: protein_name,
        7: psm_pep, 8: psm_calc_mz, 9: psm_obs_mz, 10: psm_is_decoy,
        11: psm_scan, 12: psm_id_run, 13: psm_rt,
        14: pg_global_qvalue, 15: gene_name, 16: modifications_json
        """
        records: list[dict] = []
        _enzyme = self._enzyme_name
        _count_mc = count_missed_cleavages
        _json_loads = json.loads
        _pep_map = getattr(self, "_pep_protein_map", {})

        for row in rows:
            try:
                peptidoform = str(row[0]) if row[0] else ""
                sequence = str(row[1]) if row[1] else ""
                run_file_name = str(row[2]) if row[2] else ""
                charge = int(row[3]) if row[3] is not None else 0
                intensity_val = float(row[4]) if row[4] is not None else 0.0
                rt = float(row[5]) if row[5] is not None else None
                protein_name = str(row[6]) if row[6] else ""

                # PSM info
                pep = float(row[7]) if row[7] is not None else None
                calc_mz = float(row[8]) if row[8] is not None else None
                obs_mz = float(row[9]) if row[9] is not None else None
                is_decoy = bool(row[10]) if row[10] is not None else False
                scan = [int(row[11])] if row[11] is not None else []
                id_run = str(row[12]) if row[12] else None
                psm_rt = float(row[13]) if row[13] is not None else None

                # Fall back to PSM RT
                if rt is None:
                    rt = psm_rt

                # Mass error
                mass_error_ppm = 1e6 * (obs_mz - calc_mz) / calc_mz if calc_mz and obs_mz else None

                # Protein accessions — resolve protein groups via PEP section
                acc_list = protein_name.split(";") if protein_name else []
                if len(acc_list) > 1 and sequence in _pep_map:
                    resolved = _pep_map[sequence]
                    acc_list = [resolved]
                anchor_protein = acc_list[0] if acc_list else ""

                # Protein q-value (from SQL JOIN)
                pg_qval = float(row[14]) if row[14] is not None else None

                # Gene name (from SQL JOIN)
                gene_name = str(row[15]) if row[15] else None
                gg = [gene_name] if gene_name else None

                # Modifications (from ProForma lookup)
                mods_json = row[16]
                modifications = _json_loads(mods_json) if mods_json else None

                intensities = [{"label": "LFQ", "intensity": intensity_val}]

                records.append(
                    {
                        "sequence": sequence,
                        "peptidoform": peptidoform,
                        "modifications": modifications,
                        "charge": charge,
                        "posterior_error_probability": pep,
                        "is_decoy": is_decoy,
                        "calculated_mz": calc_mz or 0.0,
                        "observed_mz": obs_mz or 0.0,
                        "mass_error_ppm": mass_error_ppm,
                        "missed_cleavages": (_count_mc(sequence, _enzyme) if _enzyme else None),
                        "additional_scores": None,
                        "predicted_rt": None,
                        "run_file_name": run_file_name,
                        "cv_params": None,
                        "scan": scan,
                        "rt": rt,
                        "ion_mobility": None,
                        "intensities": intensities,
                        "additional_intensities": None,
                        "pg_accessions": (
                            [
                                {
                                    "accession": a,
                                    "start": None,
                                    "end": None,
                                    "pre": None,
                                    "post": None,
                                }
                                for a in acc_list
                            ]
                            if acc_list
                            else None
                        ),
                        "anchor_protein": anchor_protein,
                        "unique": len(acc_list) <= 1,
                        "pg_global_qvalue": pg_qval,
                        "pg_positions": None,
                        "ion_mobility_start": None,
                        "ion_mobility_stop": None,
                        "gg_accessions": gg,
                        "gg_names": gg,
                        "id_run_file_name": id_run,
                        "rt_start": None,
                        "rt_stop": None,
                    }
                )
            except Exception as e:
                self.logger.debug(f"Skipping feature row: {e}")

        return records

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    # Standard TMT channel index → label mappings
    _TMT_CHANNEL_MAP: dict[int, str] = {
        1: "TMT126",
        2: "TMT127N",
        3: "TMT127C",
        4: "TMT128N",
        5: "TMT128C",
        6: "TMT129N",
        7: "TMT129C",
        8: "TMT130N",
        9: "TMT130C",
        10: "TMT131N",
        11: "TMT131C",
        # TMT16plex extensions
        12: "TMT132N",
        13: "TMT132C",
        14: "TMT133N",
        15: "TMT133C",
        16: "TMT134N",
        # TMT18plex extensions
        17: "TMT134C",
        18: "TMT135N",
    }

    def _detect_experiment_type(self) -> str:
        """Detect experiment type from MSstats Channel column."""
        try:
            channels = self._conn.execute('SELECT DISTINCT "Channel" FROM msstats LIMIT 20').fetchall()
            channel_vals = [str(c[0]).upper() for c in channels if c[0]]
            if any("TMT" in c for c in channel_vals):
                return "TMT"
            if any("ITRAQ" in c for c in channel_vals):
                return "iTRAQ"
            # Numeric channels with >1 distinct value indicate isobaric labeling
            if len(channel_vals) > 1:
                try:
                    int_vals = [int(float(v)) for v in channel_vals]
                    if all(1 <= v <= 18 for v in int_vals):
                        return "TMT"
                except (ValueError, TypeError):
                    pass
        except Exception:
            self.logger.debug("Failed to detect labeling type from msstats channels", exc_info=True)
        return "LFQ"

    def _build_psm_lookup(self, ms_runs: dict[int, str]) -> dict[tuple, dict]:
        """Build a lookup from (run_file_name, peptidoform, charge) -> PSM info.

        Uses a DuckDB SQL pre-aggregation to extract only the columns and
        first-occurrence rows needed, avoiding a full Python-side iteration
        over all 500K+ PSM rows.
        """
        lookup: dict[tuple, dict] = {}

        # Detect which CV column names are present for peptidoform / decoy
        actual_cols = {
            c[0]
            for c in self._conn.execute("SELECT column_name FROM information_schema.columns WHERE table_name='psms'").fetchall()
        }

        # Peptidoform column (CV_PEPTIDOFORM_SEQUENCE = MS:1000889)
        pf_lo, pf_hi = (
            "opt_global_cv_ms:1000889_peptidoform_sequence",
            "opt_global_cv_MS:1000889_peptidoform_sequence",
        )
        if pf_lo in actual_cols:
            pf_col = pf_lo
        elif pf_hi in actual_cols:
            pf_col = pf_hi
        else:
            pf_col = None

        # Decoy column (CV_DECOY_PEPTIDE = MS:1002217)
        dec_lo, dec_hi = (
            "opt_global_cv_ms:1002217_decoy_peptide",
            "opt_global_cv_MS:1002217_decoy_peptide",
        )
        if dec_lo in actual_cols:
            dec_col = dec_lo
        elif dec_hi in actual_cols:
            dec_col = dec_hi
        else:
            dec_col = None

        # PEP column
        pep_col = None
        for candidate in [
            "opt_global_posterior_error_probability_score",
            "opt_global_Posterior_Error_Probability_score",
        ]:
            if candidate in actual_cols:
                pep_col = candidate
                break

        try:
            # Build SQL to extract exactly the needed columns
            pf_expr = (
                sql_build(
                    "CAST($col AS VARCHAR)",
                    col=validate_identifier(pf_col),
                )
                if pf_col
                else "''"
            )
            dec_expr = (
                sql_build(
                    "CAST($col AS VARCHAR)",
                    col=validate_identifier(dec_col),
                )
                if dec_col
                else "'0'"
            )
            pep_expr = (
                sql_build(
                    "TRY_CAST($col AS DOUBLE)",
                    col=validate_identifier(pep_col),
                )
                if pep_col
                else "NULL"
            )
            rt_expr = "TRY_CAST(retention_time AS DOUBLE)" if "retention_time" in actual_cols else "NULL"

            stmt = sql_build(
                """
                SELECT
                    spectra_ref,
                    $pf_expr AS peptidoform,
                    CAST(charge AS VARCHAR) AS charge,
                    $pep_expr AS pep,
                    TRY_CAST(calc_mass_to_charge AS DOUBLE) AS calc_mz,
                    TRY_CAST(exp_mass_to_charge AS DOUBLE) AS obs_mz,
                    $dec_expr AS is_decoy_raw,
                    CAST(accession AS VARCHAR) AS accession,
                    $rt_expr AS rt
                FROM psms
                """,
                pf_expr=pf_expr,
                dec_expr=dec_expr,
                pep_expr=pep_expr,
                rt_expr=rt_expr,
            )
            rows = self._conn.execute(stmt).fetchall()

            for (
                spectra_ref,
                peptidoform,
                charge,
                pep,
                calc_mz,
                obs_mz,
                is_decoy_raw,
                accession,
                rt_val,
            ) in rows:
                spectra_ref = str(spectra_ref) if spectra_ref else ""
                run_file_raw = resolve_run_file(spectra_ref, ms_runs) or ""
                run_file = run_file_raw.split(".")[0] if run_file_raw else ""
                peptidoform = str(peptidoform) if peptidoform else ""
                charge = str(charge) if charge else "0"

                key = (run_file, peptidoform, charge)
                if key not in lookup:
                    is_decoy = str(is_decoy_raw).strip() == "1" if is_decoy_raw else False
                    scan = parse_scan_numbers(spectra_ref)
                    lookup[key] = {
                        "pep": float(pep) if pep is not None else None,
                        "calculated_mz": (float(calc_mz) if calc_mz is not None else None),
                        "observed_mz": float(obs_mz) if obs_mz is not None else None,
                        "is_decoy": is_decoy,
                        "accession": str(accession) if accession else "",
                        "scan": scan,
                        "id_run_file_name": run_file,
                        "rt": float(rt_val) if rt_val is not None else None,
                    }
        except Exception as e:
            self.logger.warning("Could not build PSM lookup: %s", e, exc_info=True)

        return lookup

    def _build_protein_qvalue_map(self) -> dict[str, float]:
        """Build protein accession -> global q-value."""
        try:
            # Detect the actual score column name (may contain brackets)
            cols = {
                c[0]
                for c in self._conn.execute(
                    "SELECT column_name FROM information_schema.columns WHERE table_name='proteins'"
                ).fetchall()
            }
            score_col = None
            for candidate in [
                "best_search_engine_score[1]",
                "best_search_engine_score_1",
            ]:
                if candidate in cols:
                    score_col = candidate
                    break
            if not score_col:
                return {}
            q_sc = validate_identifier(score_col)
            rows = self._conn.execute(
                sql_build(
                    """SELECT accession, $sc
                FROM proteins
                WHERE accession IS NOT NULL
                  AND $sc IS NOT NULL""",
                    sc=q_sc,
                )
            ).fetchall()
            return {str(acc): float(qval) for acc, qval in rows}
        except Exception:
            return {}

    def _build_protein_gene_map(self) -> dict[str, list[str]]:
        """Build protein accession -> gene symbol(s) from mzTab proteins table.

        Gene names are extracted from the UniProt description field using
        the ``GN=`` tag, mirroring the logic in
        ``QuantmsPgAdapter._parse_protein_row``.
        """
        gene_map: dict[str, list[str]] = {}
        try:
            rows = self._conn.execute("SELECT accession, description FROM proteins WHERE accession IS NOT NULL").fetchall()
            for acc, desc in rows:
                acc_str = str(acc).strip()
                if not acc_str or acc_str == "null":
                    continue
                if desc and str(desc) != "null":
                    gn = re.search(r"GN=([^\s]+)", str(desc))
                    if gn:
                        gene_map[acc_str] = [gn.group(1)]
        except Exception:
            self.logger.debug("Failed to extract gene map from proteins table", exc_info=True)
        return gene_map

    def _iter_feature_batches(self, file_batch_size: int):
        """Yield DataFrames of MSstats data grouped by reference files."""
        try:
            ref_col = self._detect_ref_column()
            q_ref = validate_identifier(ref_col)
            refs = self._conn.execute(
                sql_build(
                    "SELECT DISTINCT $col FROM msstats ORDER BY $col",
                    col=q_ref,
                )
            ).fetchall()
            ref_list = [r[0] for r in refs]

            for i in range(0, len(ref_list), file_batch_size):
                batch_refs = ref_list[i : i + file_batch_size]
                placeholders = ", ".join(["?" for _ in batch_refs])
                df = self._conn.execute(
                    sql_build(
                        "SELECT * FROM msstats WHERE $col IN ($ph)",
                        col=q_ref,
                        ph=placeholders,
                    ),
                    batch_refs,
                ).df()
                if not df.empty:
                    yield df
        except Exception as e:
            self.logger.error(f"Error iterating feature batches: {e}")

    def _detect_ref_column(self) -> str:
        """Detect the reference/run column name in MSstats table.

        Uses FIELD_MAPPINGS candidates if available, falls back to common names.
        """
        cols = self._conn.execute("SELECT column_name FROM information_schema.columns WHERE table_name='msstats'").fetchall()
        col_names = {c[0] for c in cols}
        # Try FIELD_MAPPINGS candidates first
        for candidate in _FEATURE_MAP.get("run_file_name", []):
            if candidate in col_names:
                return candidate
        # Fallback
        for candidate in ["Reference", "reference", "Run", "run"]:
            if candidate in col_names:
                return candidate
        return "Reference"

    def _transform_batch(
        self,
        df: pd.DataFrame,
        psm_lookup: dict,
        protein_qvalue_map: dict,
        protein_gene_map: dict,
        experiment_type: str,
        ms_runs: dict,
    ) -> list[dict]:
        """Transform an MSstats batch DataFrame into QPX feature records.

        Dispatches to a fast LFQ path (direct row iteration, no groupby)
        or an isobaric path (groupby to aggregate channels per feature).
        """
        col_map = resolve_columns(_FEATURE_MAP, set(df.columns))
        charge_col = col_map.get("charge", "Charge")
        if charge_col not in df.columns:
            df[charge_col] = 0

        if experiment_type == "LFQ":
            return self._transform_batch_lfq(df, col_map, psm_lookup, protein_qvalue_map, protein_gene_map)
        return self._transform_batch_isobaric(
            df,
            col_map,
            psm_lookup,
            protein_qvalue_map,
            protein_gene_map,
            experiment_type,
        )

    # ------ LFQ fast path (1 row per feature, no groupby) ------

    def _transform_batch_lfq(
        self,
        df: pd.DataFrame,
        col_map: dict,
        psm_lookup: dict,
        protein_qvalue_map: dict,
        protein_gene_map: dict,
    ) -> list[dict]:
        """Build feature records for LFQ data by iterating rows directly.

        For LFQ, each MSstats row is a unique (peptidoform, protein, run, charge)
        combination, so groupby is pure overhead.  We use numpy column arrays
        for fast element access and avoid ``to_dict("records")`` entirely.
        """
        records: list[dict] = []

        pf_col = col_map.get("peptidoform", "PeptideSequence")
        prot_col = col_map.get("pg_accessions", "ProteinName")
        ref_col = col_map.get("run_file_name", "Reference")
        charge_col = col_map.get("charge", "Charge")
        intensity_col = col_map.get("intensity", "Intensity")
        rt_col = col_map.get("rt", "RetentionTime")

        # Pre-extract numpy arrays for fast element access
        pf_arr = df[pf_col].values
        prot_arr = df[prot_col].values
        ref_arr = df[ref_col].values
        charge_arr = df[charge_col].values
        int_arr = df[intensity_col].values
        rt_arr = df[rt_col].values if rt_col in df.columns else None

        _sub = re.sub
        _from_proforma = from_proforma
        _safe_float = safe_float
        mods_meta = self._modifications_meta
        _proforma_cache: dict[tuple[str, str], list | None] = {}

        n = len(df)
        for i in range(n):
            try:
                peptidoform = str(pf_arr[i]) if pf_arr[i] is not None and pd.notna(pf_arr[i]) else ""
                protein_name = str(prot_arr[i]) if prot_arr[i] is not None and pd.notna(prot_arr[i]) else ""
                run_file_name = str(ref_arr[i]).split(".")[0] if ref_arr[i] is not None and pd.notna(ref_arr[i]) else ""

                charge_raw = charge_arr[i]
                charge = int(float(charge_raw)) if charge_raw is not None and charge_raw != "" and pd.notna(charge_raw) else 0

                sequence = _sub(r"[^A-Z]", "", peptidoform.upper()) if peptidoform else ""

                if peptidoform and peptidoform != sequence:
                    _cache_key = (peptidoform, sequence)
                    if _cache_key in _proforma_cache:
                        peptidoform, modifications = _proforma_cache[_cache_key]
                    else:
                        peptidoform, modifications = _from_proforma(
                            peptidoform,
                            sequence,
                            meta=mods_meta,
                        )
                        _proforma_cache[_cache_key] = (peptidoform, modifications)
                else:
                    modifications = None

                intensity_val = _safe_float(int_arr[i]) or 0.0
                intensities = [{"label": "LFQ", "intensity": float(intensity_val)}]

                rt = _safe_float(rt_arr[i]) if rt_arr is not None else None

                psm_key = (run_file_name, peptidoform, str(charge))
                psm_info = psm_lookup.get(psm_key, {})

                # Fall back to PSM RT when MSstats has no RT
                if rt is None:
                    rt = psm_info.get("rt")

                _calc = psm_info.get("calculated_mz")
                _obs = psm_info.get("observed_mz")
                mass_error_ppm = 1e6 * (_obs - _calc) / _calc if _calc and _obs else None

                acc_list = protein_name.split(";") if protein_name else []
                anchor_protein = acc_list[0] if acc_list else ""

                records.append(
                    {
                        "sequence": sequence,
                        "peptidoform": peptidoform,
                        "modifications": modifications,
                        "charge": charge,
                        "posterior_error_probability": psm_info.get("pep"),
                        "is_decoy": psm_info.get("is_decoy", False),
                        "calculated_mz": _calc or 0.0,
                        "observed_mz": _obs or 0.0,
                        "mass_error_ppm": mass_error_ppm,
                        "missed_cleavages": (count_missed_cleavages(sequence, self._enzyme_name) if self._enzyme_name else None),
                        "additional_scores": None,
                        "predicted_rt": None,
                        "run_file_name": run_file_name,
                        "cv_params": None,
                        "scan": psm_info.get("scan", []),
                        "rt": rt,
                        "ion_mobility": None,
                        "intensities": intensities,
                        "additional_intensities": None,
                        "pg_accessions": (
                            [
                                {
                                    "accession": a,
                                    "start": None,
                                    "end": None,
                                    "pre": None,
                                    "post": None,
                                }
                                for a in acc_list
                            ]
                            if acc_list
                            else None
                        ),
                        "anchor_protein": anchor_protein,
                        "unique": len(acc_list) <= 1,
                        "pg_global_qvalue": protein_qvalue_map.get(anchor_protein),
                        "pg_positions": None,
                        "ion_mobility_start": None,
                        "ion_mobility_stop": None,
                        "gg_accessions": protein_gene_map.get(anchor_protein),
                        "gg_names": protein_gene_map.get(anchor_protein),
                        "id_run_file_name": psm_info.get("id_run_file_name"),
                        "rt_start": None,
                        "rt_stop": None,
                    }
                )
            except Exception as e:
                self.logger.debug(f"Skipping feature row {i}: {e}")

        return records

    # ------ Isobaric path (TMT / iTRAQ -- groupby to aggregate channels) ------

    def _transform_batch_isobaric(
        self,
        df: pd.DataFrame,
        col_map: dict,
        psm_lookup: dict,
        protein_qvalue_map: dict,
        protein_gene_map: dict,
        experiment_type: str,
    ) -> list[dict]:
        """Build feature records for isobaric (TMT/iTRAQ) data.

        Uses groupby to aggregate channel intensities per feature, but
        accesses column values via ``.values`` arrays instead of
        ``to_dict("records")``.
        """
        records: list[dict] = []

        pf_col = col_map.get("peptidoform", "PeptideSequence")
        prot_col = col_map.get("pg_accessions", "ProteinName")
        ref_col = col_map.get("run_file_name", "Reference")
        charge_col = col_map.get("charge", "Charge")
        intensity_col = col_map.get("intensity", "Intensity")
        channel_col = col_map.get("channel", "Channel")
        rt_col = col_map.get("rt", "RetentionTime")

        grouping = [pf_col, prot_col, ref_col, charge_col]

        _sub = re.sub
        _from_proforma = from_proforma
        _safe_float = safe_float
        _tmt_map = self._TMT_CHANNEL_MAP
        mods_meta = self._modifications_meta
        is_tmt = experiment_type == "TMT"
        _proforma_cache: dict[tuple[str, str], list | None] = {}

        for group_key, group_data in df.groupby(grouping, dropna=False):
            try:
                peptidoform, protein_name, run_file_name, charge = group_key

                peptidoform = str(peptidoform) if peptidoform else ""
                protein_name = str(protein_name) if protein_name else ""
                run_file_name = str(run_file_name).split(".")[0] if run_file_name else ""
                charge = int(float(charge)) if charge not in (None, "", "null") else 0

                sequence = _sub(r"[^A-Z]", "", peptidoform.upper()) if peptidoform else ""

                if peptidoform and peptidoform != sequence:
                    _cache_key = (peptidoform, sequence)
                    if _cache_key in _proforma_cache:
                        peptidoform, modifications = _proforma_cache[_cache_key]
                    else:
                        peptidoform, modifications = _from_proforma(
                            peptidoform,
                            sequence,
                            meta=mods_meta,
                        )
                        _proforma_cache[_cache_key] = (peptidoform, modifications)
                else:
                    modifications = None

                # Build intensities from channel values (no to_dict)
                ch_vals = group_data[channel_col].values
                int_vals = group_data[intensity_col].values
                intensities = []
                for j in range(len(ch_vals)):
                    ch_raw = ch_vals[j]
                    if is_tmt and ch_raw is not None and pd.notna(ch_raw) and ch_raw != "":
                        try:
                            label = _tmt_map.get(int(float(ch_raw)), str(ch_raw))
                        except (ValueError, TypeError):
                            label = str(ch_raw)
                    else:
                        label = str(ch_raw) if ch_raw is not None and pd.notna(ch_raw) and ch_raw != "" else "LFQ"
                    iv = _safe_float(int_vals[j]) or 0.0
                    intensities.append({"label": label, "intensity": float(iv)})

                rt = _safe_float(group_data[rt_col].values[0]) if rt_col in group_data.columns else None

                psm_key = (run_file_name, peptidoform, str(charge))
                psm_info = psm_lookup.get(psm_key, {})

                # Fall back to PSM RT when MSstats has no RT
                if rt is None:
                    rt = psm_info.get("rt")

                _calc = psm_info.get("calculated_mz")
                _obs = psm_info.get("observed_mz")
                mass_error_ppm = 1e6 * (_obs - _calc) / _calc if _calc and _obs else None

                acc_list = protein_name.split(";") if protein_name else []
                anchor_protein = acc_list[0] if acc_list else ""

                records.append(
                    {
                        "sequence": sequence,
                        "peptidoform": peptidoform,
                        "modifications": modifications,
                        "charge": charge,
                        "posterior_error_probability": psm_info.get("pep"),
                        "is_decoy": psm_info.get("is_decoy", False),
                        "calculated_mz": _calc or 0.0,
                        "observed_mz": _obs or 0.0,
                        "mass_error_ppm": mass_error_ppm,
                        "missed_cleavages": (count_missed_cleavages(sequence, self._enzyme_name) if self._enzyme_name else None),
                        "additional_scores": None,
                        "predicted_rt": None,
                        "run_file_name": run_file_name,
                        "cv_params": None,
                        "scan": psm_info.get("scan", []),
                        "rt": rt,
                        "ion_mobility": None,
                        "intensities": intensities or None,
                        "additional_intensities": None,
                        "pg_accessions": (
                            [
                                {
                                    "accession": a,
                                    "start": None,
                                    "end": None,
                                    "pre": None,
                                    "post": None,
                                }
                                for a in acc_list
                            ]
                            if acc_list
                            else None
                        ),
                        "anchor_protein": anchor_protein,
                        "unique": len(acc_list) <= 1,
                        "pg_global_qvalue": protein_qvalue_map.get(anchor_protein),
                        "pg_positions": None,
                        "ion_mobility_start": None,
                        "ion_mobility_stop": None,
                        "gg_accessions": protein_gene_map.get(anchor_protein),
                        "gg_names": protein_gene_map.get(anchor_protein),
                        "id_run_file_name": psm_info.get("id_run_file_name"),
                        "rt_start": None,
                        "rt_stop": None,
                    }
                )
            except Exception as e:
                self.logger.debug(f"Skipping feature group: {e}")

        return records

    @staticmethod
    def _detect_msstats_columns(df: pd.DataFrame) -> dict[str, str]:
        """Detect the actual MSstats column names in a DataFrame.

        .. deprecated:: Use ``resolve_columns(FIELD_MAPPINGS["feature"], ...)`` instead.
           Kept for backward compatibility with external callers.
        """
        return resolve_columns(_FEATURE_MAP, set(df.columns))
