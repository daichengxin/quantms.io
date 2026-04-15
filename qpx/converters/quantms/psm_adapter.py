"""QuantMS PSM adapter -- mzTab PSM section to psm.parquet.

Loads an mzTab file into DuckDB, transforms the PSM section into the
new ``PsmSchema`` using SQL, and streams the results through
``PsmWriter``.

Key schema changes handled here:
    - ``reference_file_name`` -> ``run_file_name``
    - ``precursor_charge`` (int32) -> ``charge`` (int16)
    - ``is_decoy`` (int32, 0/1) -> ``is_decoy`` (bool)
    - ``scan`` (string) -> ``scan`` (list<int32>)
"""

from __future__ import annotations

import logging
import math
from typing import Optional

from qpx.converters.base import BaseConverter, resolve_columns
from qpx.converters.mappings import get_extra, get_field_mappings
from qpx.converters.mztab import (
    extract_modifications,
    extract_ms_runs,
    extract_score_names,
    load_mztab_sections,
)
from qpx.converters.ptm import from_proforma
from qpx.converters.utils import (
    get_cv_value,
    parse_scan_numbers,
    resolve_run_file,
    safe_float,
)
from qpx.core.cleavage import count_missed_cleavages
from qpx.core.cv_terms import CV_DECOY_PEPTIDE, CV_PEPTIDOFORM_SEQUENCE
from qpx.core.scores import is_higher_better, normalize_score_name
from qpx.core.sql import sql_build
from qpx.writers.psm import PsmWriter

logger = logging.getLogger(__name__)

# Derive field map from central mappings
_PSM_MAP = get_field_mappings("quantms", "psm")
_PHOSPHO = get_extra("quantms", "phospho_site_columns")


def _parse_site_probability_string(raw: str, score_name: str) -> dict[int, list[dict]]:
    """Parse per-position site probability strings like ``'S3:0.95;T5:0.03'``.

    Returns a dict mapping position (1-indexed) to a list of score dicts,
    or an empty dict if *raw* does not look like a per-position string.
    """
    result: dict[int, list[dict]] = {}
    for part in raw.split(";"):
        part = part.strip()
        if ":" not in part:
            continue
        pos_str, val_str = part.split(":", 1)
        # Extract numeric position (e.g., "S3" -> 3)
        digits = "".join(c for c in pos_str if c.isdigit())
        if not digits:
            continue
        pos = int(digits)
        try:
            val = float(val_str)
        except ValueError:
            continue
        result.setdefault(pos, []).append({"score_name": score_name, "score_value": val, "higher_better": True})
    return result


class QuantmsPsmAdapter(BaseConverter):
    """Convert QuantMS mzTab PSM data to ``psm.parquet``.

    Usage::

        with QuantmsPsmAdapter() as adapter:
            adapter.convert(
                mztab_path="results.mzTab",
                output_path="psm.parquet",
            )
    """

    def convert(
        self,
        mztab_path: str,
        output_path: str,
        chunksize: int = 500_000,
        creator: str = "quantms",
        enzyme_name: str | None = None,
    ) -> None:
        """Run the mzTab PSM -> psm.parquet conversion.

        Args:
            mztab_path: Path to the mzTab file.
            output_path: Destination path for the Parquet output.
            chunksize: Rows per Arrow batch when streaming.
            creator: Creator tag in Parquet footer metadata.
            enzyme_name: Enzyme name for computing missed cleavages
                when the mzTab file does not report them.
        """
        self._enzyme_name = enzyme_name
        # Step 1: Load mzTab sections into DuckDB (skip if already loaded)
        if not self._table_exists("psms"):
            load_mztab_sections(self._conn, mztab_path)

        # Step 2: Resolve column mappings against actual PSM table columns
        actual_cols = {
            c[0]
            for c in self._conn.execute("SELECT column_name FROM information_schema.columns WHERE table_name='psms'").fetchall()
        }
        self._resolved = resolve_columns(_PSM_MAP, actual_cols)

        # Detect phospho site localization columns present in PSM table
        self._phospho_cols = {col: score_name for col, score_name in (_PHOSPHO or {}).items() if col in actual_cols}
        if self._phospho_cols:
            self.logger.info(f"Detected phospho site columns: {list(self._phospho_cols.keys())}")

        # Step 3: Extract auxiliary lookups
        ms_runs = extract_ms_runs(self._conn)
        modifications_meta = extract_modifications(self._conn)
        score_names = extract_score_names(self._conn)

        # Pre-compute score metadata (name normalization + direction) once
        # instead of per-row, avoiding repeated OBO queries.
        self._psm_score_meta: dict[int, tuple[str, bool]] = {}
        for idx, raw_name in score_names.get("psms", {}).items():
            normalized = normalize_score_name(raw_name)
            higher = is_higher_better(normalized)
            self._psm_score_meta[idx] = (normalized, higher)

        # Step 4: Build protein q-value map
        protein_qvalue_map = self._build_protein_qvalue_map()

        # Step 5: Stream PSM rows and transform
        self.logger.info("Transforming PSM rows ...")

        with PsmWriter(output_path, creator=creator, compression=self._compression) as writer:
            for batch in self._iter_psm_batches(chunksize):
                records = self._transform_batch(
                    batch,
                    ms_runs,
                    modifications_meta,
                    score_names,
                    protein_qvalue_map,
                )
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)

        self.logger.info(f"PSM conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _build_protein_qvalue_map(self) -> dict[str, float]:
        """Build accession -> global q-value mapping from protein table."""
        try:
            rows = self._conn.execute("""
                SELECT accession, best_search_engine_score_1
                FROM proteins
                WHERE accession IS NOT NULL
                  AND best_search_engine_score_1 IS NOT NULL
                """).fetchall()
            return {acc: float(qval) for acc, qval in rows}
        except Exception:
            return {}

    def _iter_psm_batches(self, chunksize: int):
        """Yield batches from the ``psms`` table."""
        total = self._conn.execute("SELECT COUNT(*) FROM psms").fetchone()[0]
        offset = 0
        while offset < total:
            df = self._conn.execute(
                sql_build(
                    "SELECT * FROM psms LIMIT $lim OFFSET $off",
                    lim=str(int(chunksize)),
                    off=str(int(offset)),
                )
            ).df()
            if df.empty:
                break
            yield df
            offset += chunksize

    def _transform_batch(
        self,
        df,
        ms_runs: dict[int, str],
        modifications_meta: dict,
        score_names: dict,
        protein_qvalue_map: dict[str, float],
    ) -> list[dict]:
        """Transform a pandas DataFrame of raw mzTab PSM rows into QPX records."""
        records: list[dict] = []
        skipped = 0

        for row in df.to_dict("records"):
            try:
                rec = self._transform_row(row, ms_runs, modifications_meta, score_names, protein_qvalue_map)
                if rec:
                    records.append(rec)
            except Exception as e:
                skipped += 1
                self.logger.debug(f"Skipping PSM row: {e}")

        if skipped:
            total = skipped + len(records)
            self.logger.warning(
                "Skipped %d / %d rows (%.1f%%) in batch",
                skipped,
                total,
                100 * skipped / total if total else 0,
            )

        return records

    def _transform_row(
        self,
        row,
        ms_runs: dict[int, str],
        modifications_meta: dict,
        score_names: dict,
        protein_qvalue_map: dict[str, float],
    ) -> Optional[dict]:
        """Transform a single PSM row dict into QPX schema."""

        r = self._resolved  # shorthand for resolved column mappings

        spectra_ref = str(row.get("spectra_ref", ""))

        # --- Core identification ---
        sequence = str(row.get(r.get("sequence", "sequence"), ""))
        # mzTab column embeds CV_PEPTIDOFORM_SEQUENCE (MS:1000889)
        peptidoform_raw = get_cv_value(row, CV_PEPTIDOFORM_SEQUENCE, "peptidoform_sequence", "")
        peptidoform = str(peptidoform_raw) if peptidoform_raw else sequence

        charge_raw = row.get(r.get("charge", "charge"), 0)
        charge = int(float(charge_raw)) if charge_raw not in (None, "", "null") else 0

        # --- Decoy flag (bool) ---
        # mzTab column embeds CV_DECOY_PEPTIDE (MS:1002217)
        is_decoy_raw = get_cv_value(row, CV_DECOY_PEPTIDE, "decoy_peptide", "0")
        is_decoy = str(is_decoy_raw).strip() == "1"

        # --- Scan (list<int32>) ---
        scan = parse_scan_numbers(spectra_ref)

        # --- Run file name ---
        run_file_name = resolve_run_file(spectra_ref, ms_runs) or ""

        # --- m/z values ---
        observed_mz = safe_float(row.get(r.get("observed_mz", "exp_mass_to_charge")))
        calculated_mz = safe_float(row.get(r.get("calculated_mz", "calc_mass_to_charge")))

        # --- Mass error (ppm) ---
        if observed_mz and calculated_mz:
            mass_error_ppm = 1e6 * (observed_mz - calculated_mz) / calculated_mz
        else:
            mass_error_ppm = None

        # --- Retention time ---
        rt = safe_float(row.get(r.get("rt", "retention_time")))

        # --- PEP ---
        pep_col = r.get(
            "posterior_error_probability",
            "opt_global_Posterior_Error_Probability_score",
        )
        pep = safe_float(
            row.get(
                "opt_global_posterior_error_probability_score",
                row.get(pep_col),
            )
        )

        # --- Predicted RT ---
        predicted_rt = safe_float(row.get("opt_global_predicted_rt"))

        # --- Protein accessions ---
        accession_raw = str(row.get("accession", ""))
        protein_accessions = sorted(accession_raw.split(",")) if accession_raw else []

        # --- Additional scores ---
        additional_scores = self._build_additional_scores(row, score_names)

        # Append protein global q-value if available
        # Use the minimum q-value among individual protein accessions
        pg_qvals = [protein_qvalue_map[acc] for acc in protein_accessions if acc in protein_qvalue_map]
        pg_qval = min(pg_qvals) if pg_qvals else None
        if pg_qval is not None:
            additional_scores.append(
                {
                    "score_name": "protein_global_qvalue",
                    "score_value": float(pg_qval),
                    "higher_better": False,
                }
            )

        # global_qvalue from PSM table
        global_qvalue = safe_float(row.get("opt_global_q-value", row.get("global_qvalue")))
        if global_qvalue is not None:
            additional_scores.append(
                {
                    "score_name": "global_qvalue",
                    "score_value": global_qvalue,
                    "higher_better": False,
                }
            )

        # --- Flanking AAs and peptide positions (standard mzTab fields) ---
        pre = row.get("pre")
        post = row.get("post")
        start = row.get("start")
        end = row.get("end")

        # --- CV params ---
        cv_params = []

        def _is_valid(val):
            if val in (None, "", "null"):
                return False
            if isinstance(val, float) and math.isnan(val):
                return False
            return True

        if _is_valid(pre):
            cv_params.append({"cv_name": "pre", "cv_value": str(pre)})
        if _is_valid(post):
            cv_params.append({"cv_name": "post", "cv_value": str(post)})
        if _is_valid(start):
            cv_params.append({"cv_name": "start", "cv_value": str(start)})
        if _is_valid(end):
            cv_params.append({"cv_name": "end", "cv_value": str(end)})

        consensus = row.get(
            "opt_global_consensus_support",
            row.get("consensus_support"),
        )
        if consensus not in (None, "", "null"):
            cv_params.append({"cv_name": "consensus_support", "cv_value": str(consensus)})

        # --- Phospho site localization scores ---
        # Scan for detected phospho columns; attach as per-PSM additional_scores
        # and build site_scores for modification positions
        site_scores: dict[int, list[dict]] | None = None
        for col, score_name in self._phospho_cols.items():
            raw_val = row.get(col)
            if raw_val in (None, "", "null"):
                continue
            raw_str = str(raw_val)
            # Try parsing as per-position string (e.g. "S3:0.95;T5:0.03")
            parsed_positions = _parse_site_probability_string(raw_str, score_name)
            if parsed_positions:
                if site_scores is None:
                    site_scores = {}
                for pos, scores_list in parsed_positions.items():
                    site_scores.setdefault(pos, []).extend(scores_list)
            # Also add as a PSM-level score if it's a single float
            val = safe_float(raw_val)
            if val is not None:
                additional_scores.append(
                    {
                        "score_name": score_name,
                        "score_value": val,
                        "higher_better": True,
                    }
                )

        # --- Modifications (structured) ---
        modifications = from_proforma(
            peptidoform,
            sequence,
            meta=modifications_meta,
            site_scores=site_scores,
        )

        return {
            "sequence": sequence,
            "peptidoform": peptidoform,
            "modifications": modifications,
            "charge": charge,
            "posterior_error_probability": pep,
            "is_decoy": is_decoy,
            "calculated_mz": calculated_mz or 0.0,
            "observed_mz": observed_mz or 0.0,
            "mass_error_ppm": mass_error_ppm,
            "additional_scores": additional_scores or None,
            "predicted_rt": predicted_rt,
            "run_file_name": run_file_name,
            "cv_params": cv_params or None,
            "scan": scan,
            "rt": rt,
            "ion_mobility": None,
            "missed_cleavages": (count_missed_cleavages(sequence, self._enzyme_name) if self._enzyme_name else None),
            "protein_accessions": protein_accessions or None,
            "mz_array": None,
            "intensity_array": None,
            "charge_array": None,
            "ion_type_array": None,
            "ion_mobility_array": None,
        }

    def _build_additional_scores(self, row, score_names: dict) -> list[dict]:
        """Build structured additional_scores from search engine score columns.

        Uses ``self._psm_score_meta`` (pre-computed in ``convert()``) to
        avoid per-row score normalization and OBO lookups.
        """
        scores: list[dict] = []
        for idx, (normalized, higher) in self._psm_score_meta.items():
            col = f"search_engine_score[{idx}]"
            # Try lower-case variant
            val = row.get(col, row.get(col.lower()))
            if val not in (None, "", "null"):
                try:
                    scores.append(
                        {
                            "score_name": normalized,
                            "score_value": float(val),
                            "higher_better": higher,
                        }
                    )
                except (ValueError, TypeError):
                    pass
        return scores
