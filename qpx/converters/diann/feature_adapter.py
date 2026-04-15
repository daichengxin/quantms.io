"""DIA-NN Feature adapter -- report to feature.parquet.

Performance strategy:
    - DuckDB VIEW over parquet for memory efficiency (no in-memory table)
    - Precursor lookup (peptidoform, mz, missed_cleavages) via DuckDB temp table JOIN
    - All scalar and nested fields (scores, intensities, cv_params) constructed in SQL
    - Only modifications + pg_accessions built in Python (requires PTM parsing)
    - Arrow tables written directly via write_table() — no Python dict intermediaries
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import pandas as pd
import pyarrow as pa

from qpx.converters.base import resolve_columns
from qpx.converters.diann.base_adapter import DiaNNBaseAdapter
from qpx.converters.diann.constants import to_modifications, to_proforma
from qpx.converters.mappings import get_field_mappings
from qpx.converters.ptm import compute_precursor_mz
from qpx.core.cleavage import count_missed_cleavages
from qpx.core.sql import sql_build, validate_identifier
from qpx.writers.feature import FeatureWriter

logger = logging.getLogger(__name__)

# Score mappings: (DIA-NN column, output score name, higher_better)
_SCORE_MAPPINGS = [
    ("Q.Value", "qvalue", False),
    ("PG.Q.Value", "pg_qvalue", False),
    ("Global.Q.Value", "global_qvalue", False),
    ("CScore", "diann_cscore", True),
    ("Evidence", "diann_evidence", True),
    ("Spectrum.Similarity", "diann_spectrum_similarity", True),
    ("Ms1.Profile.Corr", "diann_ms1_profile_corr", True),
    ("Mass.Evidence", "diann_mass_evidence", True),
    ("Averagine", "diann_averagine", True),
    ("Quantity.Quality", "diann_quantity_quality", True),
    ("Lib.Q.Value", "diann_lib_qvalue", False),
    ("Lib.PG.Q.Value", "diann_lib_pg_qvalue", False),
    ("Protein.Q.Value", "diann_protein_qvalue", False),
    ("Decoy.CScore", "diann_decoy_cscore", True),
    ("Decoy.Evidence", "diann_decoy_evidence", True),
    ("Translated.Quality", "diann_translated_quality", True),
    ("Translated.Q.Value", "diann_translated_qvalue", False),
]

# CV param mappings: (DIA-NN column, output cv_name)
_CV_MAPPINGS = [
    ("Proteotypic", "proteotypic"),
    ("iRT", "irt"),
    ("Predicted.iRT", "predicted_irt"),
    ("Predicted.IM", "predicted_im"),
    ("Predicted.iIM", "predicted_iim"),
    ("iIM", "iim"),
]


class DiannFeatureAdapter(DiaNNBaseAdapter):
    """Convert DIA-NN report to ``feature.parquet``.

    Usage::

        with DiannFeatureAdapter() as adapter:
            adapter.convert(
                diann_report="report.tsv",
                output_path="feature.parquet",
                mzml_info_folder="ms_info/",
                sdrf_path="sdrf.tsv",
            )
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._mz_cache: dict[tuple[str, int], float | None] = {}

    def convert(
        self,
        diann_report: str,
        output_path: str,
        mzml_info_folder: Optional[str] = None,
        sdrf_path: Optional[str] = None,
        qvalue_threshold: float = 0.01,
        file_num: int = 100,  # VIEW-based lazy IO reduces per-batch memory
        creator: str = "diann",
    ) -> None:
        """Run the DIA-NN report -> feature.parquet conversion."""
        # 1. Create VIEW over report (lazy — no data loaded into memory)
        self._load_diann_report(diann_report)

        # 2. SDRF for sample mapping and enzyme info
        enzyme_name: str | None = None
        if sdrf_path:
            self._load_sdrf_sample_map(sdrf_path)
            enzyme_name = self._load_sdrf_enzyme(sdrf_path)

        # 3. Detect report columns (once, reused for all batches)
        report_cols = self._get_report_columns()
        self._resolved = resolve_columns(get_field_mappings("diann", "feature"), report_cols)
        decoy_col = self._resolved.get("decoy")
        has_decoy_col = decoy_col is not None

        # 4. Pre-compute precursors: DuckDB temp table + Python modifications dict
        mods_lookup = self._register_precursor_lookup(enzyme_name)

        # 5. Get target Arrow schema
        target_schema = FeatureWriter._schema_class.get_arrow_schema()
        mods_type = target_schema.field("modifications").type
        pg_acc_type = target_schema.field("pg_accessions").type

        # 6. Track discovered scores upfront (all score columns present in report)
        for diann_col, score_name, _ in _SCORE_MAPPINGS:
            if diann_col in report_cols:
                self._discovered_scores.add(score_name)

        # 7. Build SQL template (once, reused for all batches)
        sql_template = self._build_batch_sql(report_cols, has_decoy_col, qvalue_threshold)

        # 8. Discover runs
        run_names = self._discover_runs(mzml_info_folder)

        # 9. Process in batches → Arrow tables → write directly
        with FeatureWriter(output_path, creator=creator, compression=self._compression) as writer:
            for i in range(0, len(run_names), file_num):
                batch_runs = run_names[i : i + file_num]
                self.logger.info(f"Processing runs {i + 1}-{min(i + file_num, len(run_names))} of {len(run_names)}")
                table = self._process_batch_arrow(
                    batch_runs,
                    sql_template,
                    mods_lookup,
                    mods_type,
                    pg_acc_type,
                    target_schema,
                    mzml_info_folder,
                )
                if table is not None:
                    writer.write_table(table)

        self.logger.info(f"DIA-NN feature conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Data loading helpers
    # ------------------------------------------------------------------

    def _get_report_columns(self) -> set[str]:
        """Return the set of column names in the DuckDB report table/view."""
        return {
            c[0]
            for c in self._conn.execute("SELECT column_name FROM information_schema.columns WHERE table_name='report'").fetchall()
        }

    def _load_sdrf_sample_map(self, sdrf_path: str) -> dict[str, str]:
        """Load SDRF and create run_file -> sample_accession mapping."""
        from qpx.core.sdrf import SDRFHandler

        handler = SDRFHandler(sdrf_path)
        return handler.get_sample_map_run()

    def _load_sdrf_enzyme(self, sdrf_path: str) -> str | None:
        """Load the first enzyme name from SDRF for missed-cleavage computation."""
        try:
            from qpx.core.sdrf import SDRFHandler

            handler = SDRFHandler(sdrf_path)
            enzymes = handler.get_enzymes()
            if enzymes:
                return str(enzymes[0])
        except Exception:
            self.logger.debug("Could not load enzyme from SDRF")
        return None

    def _discover_runs(self, mzml_info_folder: Optional[str]) -> list[str]:
        """Discover run names from ms_info files or from the report."""
        if mzml_info_folder:
            info_files = list(Path(mzml_info_folder).glob("*_ms_info.parquet"))
        else:
            info_files = []

        if info_files:
            run_names = [f.stem.replace("_ms_info", "") for f in info_files]
            self.logger.info(f"Found {len(run_names)} MS info files")
        else:
            run_col = self._resolved["run_file_name"]
            qcol = validate_identifier(run_col)
            rows = self._conn.execute(
                sql_build(
                    "SELECT DISTINCT $col FROM report ORDER BY $col",
                    col=qcol,
                )
            ).fetchall()
            run_names = [r[0].replace(".mzML", "").replace(".raw", "") for r in rows]
            self.logger.info(f"Discovered {len(run_names)} runs from report")

        return run_names

    # ------------------------------------------------------------------
    # Precursor lookup (DuckDB temp table + Python modifications dict)
    # ------------------------------------------------------------------

    def _register_precursor_lookup(self, enzyme_name: str | None) -> dict:
        """Pre-compute per-precursor values and register flat fields in DuckDB.

        Returns:
            Dict mapping (modified_seq, sequence, charge) -> modifications list.
        """
        mod_col = self._resolved["modified_sequence"]
        seq_col = self._resolved["sequence"]
        chg_col = self._resolved["charge"]

        precursors = self._conn.execute(
            sql_build(
                "SELECT DISTINCT $mc, $sc, $cc FROM report",
                mc=validate_identifier(mod_col),
                sc=validate_identifier(seq_col),
                cc=validate_identifier(chg_col),
            )
        ).fetchall()

        self.logger.info(f"Pre-computing {len(precursors):,} unique precursors")

        mods_lookup: dict[tuple[str, str, int], list | None] = {}
        rows = []

        for modified_seq, sequence, charge in precursors:
            modified_seq = str(modified_seq)
            sequence = str(sequence)
            charge = int(charge)

            peptidoform = to_proforma(modified_seq)
            modifications = to_modifications(modified_seq, sequence)

            cache_key = (modified_seq, charge)
            if cache_key not in self._mz_cache:
                self._mz_cache[cache_key] = compute_precursor_mz(modified_seq, charge)
            calculated_mz = self._mz_cache[cache_key] or 0.0

            missed = count_missed_cleavages(sequence, enzyme_name) if enzyme_name else None

            mods_lookup[(modified_seq, sequence, charge)] = modifications
            rows.append((modified_seq, sequence, charge, peptidoform, calculated_mz, missed))

        # Register flat lookup as DuckDB temp table for SQL JOIN
        self._conn.execute("DROP TABLE IF EXISTS precursor_lookup")
        self._conn.from_df(
            pd.DataFrame(
                rows,
                columns=["modified_sequence", "sequence", "charge", "peptidoform", "calculated_mz", "missed_cleavages"],
            )
        ).create("precursor_lookup")
        self.logger.info("Registered precursor_lookup temp table in DuckDB")

        return mods_lookup

    # ------------------------------------------------------------------
    # SQL template builder
    # ------------------------------------------------------------------

    def _build_batch_sql(self, report_cols: set[str], has_decoy_col: bool, qvalue_threshold: float) -> str:
        """Build the SQL query template for batch processing.

        Constructs all output columns (including nested structs) in DuckDB SQL.
        Only modifications and pg_accessions are handled in Python post-processing.
        """
        r = self._resolved
        run_col = r["run_file_name"]
        qv_col = r["qvalue"]
        pg_col = r["pg_accessions"]

        def _has(col):
            return col in report_cols

        def _sf(col):
            """Safe float: NULL/NaN → NULL, else FLOAT."""
            return f'CASE WHEN r."{col}" IS NOT NULL AND NOT isnan(CAST(r."{col}" AS DOUBLE)) THEN CAST(r."{col}" AS FLOAT) END'

        def _sd(col):
            """Safe double: NULL/NaN → NULL, else DOUBLE."""
            return f'CASE WHEN r."{col}" IS NOT NULL AND NOT isnan(CAST(r."{col}" AS DOUBLE)) THEN CAST(r."{col}" AS DOUBLE) END'

        parts = []

        # --- Flat scalar columns ---
        parts.append(f'r."{r["sequence"]}" AS sequence')
        parts.append("lk.peptidoform")
        parts.append(f'CAST(r."{r["charge"]}" AS SMALLINT) AS charge')

        # PEP
        pep_col = r.get("posterior_error_probability")
        parts.append(
            f"{_sd(pep_col)} AS posterior_error_probability" if _has(pep_col) else "NULL::DOUBLE AS posterior_error_probability"
        )

        # is_decoy
        if has_decoy_col:
            parts.append(f'CAST(r."{r["decoy"]}" AS BOOLEAN) AS is_decoy')
        else:
            parts.append(
                f"(starts_with(r.\"{pg_col}\", 'DECOY_') OR starts_with(r.\"{pg_col}\", 'decoy_') "
                f"OR starts_with(r.\"{pg_col}\", 'rev_') OR starts_with(r.\"{pg_col}\", 'REV_')) AS is_decoy"
            )

        parts.append("COALESCE(CAST(lk.calculated_mz AS FLOAT), 0.0::FLOAT) AS calculated_mz")

        # observed_mz
        omz_col = r.get("observed_mz")
        parts.append(
            f"COALESCE({_sf(omz_col)}, 0.0::FLOAT) AS observed_mz" if omz_col and _has(omz_col) else "0.0::FLOAT AS observed_mz"
        )

        # mass_error_ppm
        me_col = r.get("mass_error_ppm")
        parts.append(f"{_sf(me_col)} AS mass_error_ppm" if me_col and _has(me_col) else "NULL::FLOAT AS mass_error_ppm")

        # predicted_rt
        prt_col = r.get("predicted_rt")
        parts.append(f"{_sf(prt_col)} AS predicted_rt" if prt_col and _has(prt_col) else "NULL::FLOAT AS predicted_rt")

        # run_file_name (strip extension)
        parts.append(f"regexp_replace(r.\"{run_col}\", '\\.(mzML|raw|d)$', '') AS run_file_name")

        # scan
        ms2_col = r.get("ms2_scan")
        if ms2_col and _has(ms2_col):
            parts.append(
                f'CASE WHEN r."{ms2_col}" IS NOT NULL THEN [CAST(r."{ms2_col}" AS INTEGER)] ELSE []::INTEGER[] END AS scan'
            )
        else:
            parts.append("[]::INTEGER[] AS scan")

        # rt
        rt_col = r.get("rt")
        parts.append(f"{_sf(rt_col)} AS rt" if rt_col and _has(rt_col) else "NULL::FLOAT AS rt")

        # ion_mobility
        im_col = r.get("ion_mobility")
        parts.append(f"{_sf(im_col)} AS ion_mobility" if im_col and _has(im_col) else "NULL::FLOAT AS ion_mobility")

        # missed_cleavages
        parts.append("CAST(lk.missed_cleavages AS SMALLINT) AS missed_cleavages")

        # anchor_protein
        parts.append(f"SPLIT_PART(r.\"{pg_col}\", ';', 1) AS anchor_protein")

        # unique
        proteotypic_col = r.get("proteotypic")
        if proteotypic_col and _has(proteotypic_col):
            parts.append(f'CAST(r."{proteotypic_col}" AS INTEGER) = 1 AS "unique"')
        else:
            parts.append(f'NOT contains(r."{pg_col}", \';\') AS "unique"')

        # pg_global_qvalue
        pgqv_col = r.get("pg_global_qvalue")
        parts.append(
            f"{_sd(pgqv_col)} AS pg_global_qvalue" if pgqv_col and _has(pgqv_col) else "NULL::DOUBLE AS pg_global_qvalue"
        )

        # ion_mobility_start/stop
        for field in ("ion_mobility_start", "ion_mobility_stop"):
            col = r.get(field)
            parts.append(f"{_sf(col)} AS {field}" if col and _has(col) else f"NULL::FLOAT AS {field}")

        # gg_accessions / gg_names
        genes_col = r.get("gg_names")
        if genes_col and _has(genes_col):
            ge = f'CASE WHEN r."{genes_col}" IS NOT NULL AND r."{genes_col}" != \'\' THEN STRING_SPLIT(r."{genes_col}", \';\') ELSE NULL END'
        else:
            ge = "NULL"
        parts.append(f"{ge} AS gg_accessions")
        parts.append(f"{ge} AS gg_names")

        # id_run_file_name
        parts.append(f"regexp_replace(r.\"{run_col}\", '\\.(mzML|raw|d)$', '') AS id_run_file_name")

        # rt_start / rt_stop
        for field in ("rt_start", "rt_stop"):
            col = r.get(field)
            parts.append(f"{_sf(col)} AS {field}" if col and _has(col) else f"NULL::FLOAT AS {field}")

        # --- Nested: intensities ---
        int_col = r["intensity"]
        parts.append(f"[STRUCT_PACK(label := 'raw', intensity := COALESCE({_sf(int_col)}, 0.0::FLOAT))] AS intensities")

        # --- Nested columns built by helpers ---
        parts.append(self._build_additional_intensities_sql(r, _has))
        parts.append(self._build_additional_scores_sql(_has))
        parts.append(self._build_cv_params_sql(r, _has))

        # --- Helper columns (for Python post-processing, dropped later) ---
        mod_col = r["modified_sequence"]
        seq_col = r["sequence"]
        chg_col = r["charge"]
        parts.append(f'r."{mod_col}" AS _modified_sequence')
        parts.append(f'r."{pg_col}" AS _pg_group')

        # --- Build full query ---
        select_clause = ",\n        ".join(parts)

        sql = sql_build(
            """
        SELECT
            $select_clause
        FROM report r
        JOIN precursor_lookup lk
            ON r.$mod_col = lk.modified_sequence
            AND r.$seq_col = lk.sequence
            AND r.$chg_col = lk.charge
        WHERE r.$run_col IN ({run_placeholders})
          AND r.$qv_col < $qv_threshold
          AND r.$pg_col IS NOT NULL
        """,
            select_clause=select_clause,
            mod_col=validate_identifier(mod_col),
            seq_col=validate_identifier(seq_col),
            chg_col=validate_identifier(chg_col),
            run_col=validate_identifier(run_col),
            qv_col=validate_identifier(qv_col),
            qv_threshold=str(qvalue_threshold),
            pg_col=validate_identifier(pg_col),
        )

        return sql

    @staticmethod
    def _build_additional_intensities_sql(r: dict, _has) -> str:
        """Build SQL for the additional_intensities nested column."""
        candidates = [
            ("PG.MaxLFQ", "maxlfq"),
        ]
        if "normalize_intensity" in r:
            candidates.append((r["normalize_intensity"], "precursor_normalised"))
        if "ms1_area" in r:
            candidates.append((r["ms1_area"], "ms1_area"))
        entries, checks = [], []
        for col, name in candidates:
            if _has(col):
                chk = f'(r."{col}" IS NOT NULL AND NOT isnan(CAST(r."{col}" AS DOUBLE)))'
                checks.append(chk)
                entries.append(
                    f"CASE WHEN {chk} THEN STRUCT_PACK("
                    f"intensity_name := '{name}', "
                    f'intensity_value := CAST(r."{col}" AS FLOAT)) END'
                )
        if not entries:
            return "NULL AS additional_intensities"
        ai_list = ", ".join(entries)
        ai_check = " OR ".join(checks)
        return (
            f"CASE WHEN {ai_check} "
            f"THEN [STRUCT_PACK(label := 'raw', "
            f"intensities := LIST_FILTER([{ai_list}], x -> x IS NOT NULL))] "
            f"ELSE NULL END AS additional_intensities"
        )

    @staticmethod
    def _build_additional_scores_sql(_has) -> str:
        """Build SQL for the additional_scores nested column."""
        entries, checks = [], []
        for diann_col, score_name, higher_better in _SCORE_MAPPINGS:
            if _has(diann_col):
                hb = "true" if higher_better else "false"
                chk = f'(r."{diann_col}" IS NOT NULL AND NOT isnan(CAST(r."{diann_col}" AS DOUBLE)))'
                checks.append(chk)
                entries.append(
                    f"CASE WHEN {chk} THEN STRUCT_PACK("
                    f"score_name := '{score_name}', "
                    f'score_value := CAST(r."{diann_col}" AS DOUBLE), '
                    f"higher_better := {hb}) END"
                )
        if not entries:
            return "NULL AS additional_scores"
        sl = ", ".join(entries)
        sc = " OR ".join(checks)
        return f"CASE WHEN {sc} THEN LIST_FILTER([{sl}], x -> x IS NOT NULL) ELSE NULL END AS additional_scores"

    @staticmethod
    def _build_cv_params_sql(r: dict, _has) -> str:
        """Build SQL for the cv_params nested column."""
        entries, checks = [], []
        pqs_col = r.get("precursor_quantification_score")
        if pqs_col and _has(pqs_col):
            chk = f'(r."{pqs_col}" IS NOT NULL)'
            checks.append(chk)
            entries.append(
                f"CASE WHEN {chk} THEN STRUCT_PACK("
                f"cv_name := 'precursor_quantification_score', "
                f'cv_value := CAST(r."{pqs_col}" AS VARCHAR)) END'
            )
        for diann_col, cv_name in _CV_MAPPINGS:
            if _has(diann_col):
                chk = f'(r."{diann_col}" IS NOT NULL AND CAST(r."{diann_col}" AS VARCHAR) != \'\')'
                checks.append(chk)
                entries.append(
                    f"CASE WHEN {chk} THEN STRUCT_PACK("
                    f"cv_name := '{cv_name}', "
                    f'cv_value := CAST(r."{diann_col}" AS VARCHAR)) END'
                )
        if not entries:
            return "NULL AS cv_params"
        cl = ", ".join(entries)
        cc = " OR ".join(checks)
        return f"CASE WHEN {cc} THEN LIST_FILTER([{cl}], x -> x IS NOT NULL) ELSE NULL END AS cv_params"

    # ------------------------------------------------------------------
    # Batch processing (SQL → Arrow → write)
    # ------------------------------------------------------------------

    def _process_batch_arrow(
        self,
        runs: list[str],
        sql_template: str,
        mods_lookup: dict,
        mods_type: pa.DataType,
        pg_acc_type: pa.DataType,
        target_schema: pa.Schema,
        mzml_info_folder: Optional[str],
    ) -> Optional[pa.Table]:
        """Execute SQL batch, add Python-only columns, cast to target schema."""
        placeholders = ", ".join(["?" for _ in runs])
        sql = sql_template.replace("{run_placeholders}", placeholders)

        # Execute SQL — DuckDB pushes predicates to parquet, only reads needed data
        table = self._conn.execute(sql, runs).fetch_arrow_table()

        if len(table) == 0:
            return None

        n = len(table)

        # --- Merge scan info from mzml files if available ---
        if mzml_info_folder:
            table = self._merge_scan_info_arrow(table, mzml_info_folder)

        # --- Build modifications column from Python lookup ---
        mod_seqs = table.column("_modified_sequence").to_pylist()
        seqs = table.column("sequence").to_pylist()
        charges = table.column("charge").to_pylist()

        mods_list = [mods_lookup.get((str(ms), str(s), int(c))) for ms, s, c in zip(mod_seqs, seqs, charges)]
        mods_array = pa.array(mods_list, type=mods_type)

        # --- Build pg_accessions column ---
        pg_groups = table.column("_pg_group").to_pylist()
        pg_list = [
            [
                {"accession": acc, "start": None, "end": None, "pre": None, "post": None}
                for acc in (str(pg).split(";") if pg else [])
            ]
            or None
            for pg in pg_groups
        ]
        pg_array = pa.array(pg_list, type=pg_acc_type)

        # --- Drop helper columns ---
        table = table.drop("_modified_sequence")
        table = table.drop("_pg_group")

        # --- Assemble final table matching target schema ---
        columns = {}
        for field in target_schema:
            if field.name == "modifications":
                columns[field.name] = mods_array
            elif field.name == "pg_accessions":
                columns[field.name] = pg_array
            elif field.name in table.schema.names:
                col = table.column(field.name)
                # Cast to target type if needed (e.g. large_string → string)
                if col.type != field.type:
                    try:
                        col = col.cast(field.type, safe=False)
                    except (pa.ArrowInvalid, pa.ArrowNotImplementedError):
                        # For complex nested types, rebuild via Python round-trip
                        col = pa.array(col.to_pylist(), type=field.type)
                columns[field.name] = col
            else:
                # Missing column → NULL array
                columns[field.name] = pa.nulls(n, type=field.type)

        return pa.table(columns, schema=target_schema)

    def _merge_scan_info_arrow(self, table: pa.Table, mzml_info_folder: str) -> pa.Table:
        """Merge MS scan info from mzml files into the Arrow table."""
        df = table.to_pandas()
        merged_parts = []

        for run_name, group_df in df.groupby("run_file_name", sort=False):
            ms_info_path = next(
                Path(mzml_info_folder).glob(f"*{run_name}_ms_info.parquet"),
                None,
            )
            if ms_info_path:
                group_df = self._merge_scan_info(group_df, ms_info_path)
            merged_parts.append(group_df)

        merged_df = pd.concat(merged_parts, ignore_index=True)
        # Rebuild Arrow table preserving the original schema for non-scan columns
        return pa.Table.from_pandas(merged_df, schema=table.schema, preserve_index=False)

    def _merge_scan_info(self, run_data: pd.DataFrame, ms_info_path: Path) -> pd.DataFrame:
        """Merge MS scan info (scan, observed_mz) with report data for a single run."""
        target = pd.read_parquet(ms_info_path, columns=["rt", "scan", "precursor_mz"])
        target = target.rename(columns={"precursor_mz": "observed_mz"})
        target["rt"] = target["rt"] / 60  # Convert to minutes

        # Ensure matching types for merge_asof (SQL may produce float32, ms_info float64)
        run_data = run_data.copy()
        run_data["rt"] = run_data["rt"].astype("float64")
        target["rt"] = target["rt"].astype("float64")

        run_data = run_data.sort_values("rt")
        result = pd.merge_asof(run_data, target, on="rt", direction="nearest", suffixes=("", "_scan"))
        return result
