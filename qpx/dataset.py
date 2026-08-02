"""Dataset class — the central entry point for opening QPX datasets."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import duckdb

from qpx._version import __version__
from qpx.core.convert import QueryResult
from qpx.core.data import (
    PG,
    PSM,
    BaseStructure,
    DatasetMeta,
    Feature,
    MzSpectra,
    Ontology,
    PepMap,
    Provenance,
    Run,
    Sample,
)
from qpx.core.data.schema import ValidationIssue, ValidationResult
from qpx.core.engine import DuckDBEngine
from qpx.core.sql import escape_path, sql_build
from qpx.version import (
    QPX_SPEC_VERSION,
    QpxVersionError,
    check_pg_columns_compatible,
    check_pg_file_compatible,
)

if TYPE_CHECKING:
    import pandas as pd

_log = logging.getLogger(__name__)


class Dataset:
    """
    A QPX dataset — a directory of related Parquet/H5AD files.

    Usage:
        ds = qpx.Dataset("PXD014414/")
        ds.feature.filter("charge > 2").to_df()
        ds.pg.join(ds.run).to_df()

    Partial access:
        ds = qpx.Dataset("path/", structures=["feature", "sample", "run"])
    """

    # Data structure registry: name → (Class, file_suffix)
    _STRUCTURE_REGISTRY = {
        "psm": (PSM, ".psm.parquet"),
        "feature": (Feature, ".feature.parquet"),
        "pg": (PG, ".pg.parquet"),
        "mz": (MzSpectra, ".mz.parquet"),
        "sample": (Sample, ".sample.parquet"),
        "run": (Run, ".run.parquet"),
        "dataset": (DatasetMeta, ".dataset.parquet"),
        "ontology": (Ontology, ".ontology.parquet"),
        "provenance": (Provenance, ".provenance.parquet"),
        "pepmap": (PepMap, ".pepmap.parquet"),
    }

    def __init__(
        self,
        path: str | Path,
        structures: list[str] | None = None,
        duckdb_memory: str = "16GB",
        duckdb_threads: int = 4,
        s3_config: dict | None = None,
        file_prefix: str | None = None,
    ):
        self._is_s3 = isinstance(path, str) and path.startswith("s3://")
        self.path = path if self._is_s3 else Path(path)
        self._file_prefix = file_prefix
        self._engine = DuckDBEngine(
            memory_limit=duckdb_memory,
            threads=duckdb_threads,
            s3_config=s3_config if self._is_s3 else None,
        )

        self._structures: dict[str, BaseStructure] = {}
        self._discover_and_register(structures)

    def _discover_and_register(self, requested: list[str] | None):
        """Scan directory, find QPX files, register as DuckDB tables.

        Checks for single Parquet files first, then falls back to
        Hive-partitioned directories (e.g., feature/ with run_file_name= subdirs).
        For S3 paths, attempts to register each structure via S3 glob.
        """
        if self._is_s3:
            self._discover_s3(requested)
        else:
            self._discover_local(requested)

    def _discover_s3(self, requested: list[str] | None):
        """Register structures from S3 path."""
        for name, (cls, suffix) in self._STRUCTURE_REGISTRY.items():
            if requested and name not in requested:
                continue
            file_name = f"{self._file_prefix}{suffix}" if self._file_prefix else f"*{suffix}"
            s3_glob = f"{self.path.rstrip('/')}/{file_name}"
            try:
                # PyArrow cannot open an S3 glob to read the footer, so register
                # first and run the PG guard against the DuckDB view's columns.
                self._engine.register_s3_parquet(name, s3_glob)
                if name == "pg":
                    cols = [row[0] for row in self._engine.execute('DESCRIBE "pg"').fetchall()]
                    check_pg_columns_compatible(cols, source=s3_glob)
                self._structures[name] = cls(
                    engine=self._engine,
                    table_name=name,
                    file_path=f"{self.path}/{name}",
                )
            except QpxVersionError as exc:
                # One incompatible/old structure file must not abort the whole
                # dataset — skip it with a warning so the other structures stay
                # usable. (A direct read_pg()/PG.from_file() still raises.)
                _log.warning("Skipping incompatible structure '%s': %s", name, exc)
            except (FileNotFoundError, duckdb.IOException):
                pass  # Structure not present in S3
            except (OSError, ValueError, duckdb.Error) as exc:
                _log.warning("Failed to register '%s' from S3: %s", name, exc)

    def _discover_local(self, requested: list[str] | None):
        """Register structures from local filesystem.

        Each structure is registered independently: a bad/old/garbage file for
        one structure is skipped with a warning so the rest of the dataset stays
        usable, instead of aborting the whole Dataset construction.
        """
        for name, (cls, suffix) in self._STRUCTURE_REGISTRY.items():
            if requested and name not in requested:
                continue
            try:
                self._register_local_structure(name, cls, suffix)
            except QpxVersionError as exc:
                # One incompatible/old structure file must not abort the whole
                # dataset — skip it with a warning so the other structures stay
                # usable. (A direct read_pg()/PG.from_file() still raises.)
                _log.warning("Skipping incompatible structure '%s': %s", name, exc)
            except (OSError, ValueError, duckdb.Error) as exc:
                _log.warning("Failed to register structure '%s': %s", name, exc)

    def _register_local_structure(self, name: str, cls, suffix: str):
        """Register a single local structure (single file or Hive-partitioned dir)."""
        # Check for single file first
        pattern = f"{self._file_prefix}{suffix}" if self._file_prefix else f"*{suffix}"
        matches = sorted(self.path.glob(pattern))
        if matches:
            if len(matches) > 1:
                _log.warning(
                    "Multiple files match '%s': %s — using first: %s",
                    pattern,
                    [m.name for m in matches],
                    matches[0].name,
                )
            file_path = matches[0]  # Take first match
            if name == "pg":
                check_pg_file_compatible(file_path)
            self._engine.register_parquet(name, file_path)
            self._structures[name] = cls(
                engine=self._engine,
                table_name=name,
                file_path=file_path,
            )
        elif self._file_prefix is None:
            # Check for Hive-partitioned directory
            part_dir = self.path / name
            part_files = list(part_dir.glob("**/*.parquet")) if part_dir.is_dir() else []
            if part_files:
                if name == "pg":
                    for part_file in sorted(part_files):
                        check_pg_file_compatible(part_file)
                self._engine.register_partitioned_parquet(name, part_dir)
                self._structures[name] = cls(
                    engine=self._engine,
                    table_name=name,
                    file_path=part_dir,
                )

    # --- Data structure accessors (lazy, return None if not present) ---
    @property
    def psm(self) -> PSM | None:
        return self._structures.get("psm")

    @property
    def feature(self) -> Feature | None:
        return self._structures.get("feature")

    @property
    def pg(self) -> PG | None:
        return self._structures.get("pg")

    @property
    def mz(self) -> MzSpectra | None:
        return self._structures.get("mz")

    @property
    def sample(self) -> Sample | None:
        return self._structures.get("sample")

    @property
    def run(self) -> Run | None:
        return self._structures.get("run")

    @property
    def dataset_meta(self) -> DatasetMeta | None:
        return self._structures.get("dataset")

    @property
    def ontology(self) -> Ontology | None:
        return self._structures.get("ontology")

    @property
    def provenance(self) -> Provenance | None:
        return self._structures.get("provenance")

    @property
    def pepmap(self) -> PepMap | None:
        return self._structures.get("pepmap")

    # --- View accessors (cached) ---
    @property
    def protein_view(self):
        if not hasattr(self, "_protein_view"):
            from qpx.views.api import ProteinView

            self._protein_view = ProteinView(self)
        return self._protein_view

    @property
    def peptide_view(self):
        if not hasattr(self, "_peptide_view"):
            from qpx.views.api import PeptideView

            self._peptide_view = PeptideView(self)
        return self._peptide_view

    @property
    def identification_summary(self):
        if not hasattr(self, "_identification_summary"):
            from qpx.views.api import IdentificationSummaryView

            self._identification_summary = IdentificationSummaryView(self)
        return self._identification_summary

    @property
    def run_summary(self):
        if not hasattr(self, "_run_summary"):
            from qpx.views.api import RunSummaryView

            self._run_summary = RunSummaryView(self)
        return self._run_summary

    @property
    def modification_view(self):
        if not hasattr(self, "_modification_view"):
            from qpx.views.api import ModificationView

            self._modification_view = ModificationView(self)
        return self._modification_view

    @property
    def qc_view(self):
        if not hasattr(self, "_qc_view"):
            from qpx.views.api import QualityControlView

            self._qc_view = QualityControlView(self)
        return self._qc_view

    @property
    def sample_summary(self):
        if not hasattr(self, "_sample_summary"):
            from qpx.views.api import SampleSummaryView

            self._sample_summary = SampleSummaryView(self)
        return self._sample_summary

    @property
    def ae_view(self):
        if not hasattr(self, "_ae_view"):
            from qpx.views.api import AbsoluteExpressionView

            self._ae_view = AbsoluteExpressionView(self)
        return self._ae_view

    # --- Cross-structure queries ---
    def sql(self, query: str) -> QueryResult:
        """Execute arbitrary SQL across registered structures."""
        return QueryResult(self._engine.execute(query))

    @property
    def available_structures(self) -> list[str]:
        return list(self._structures.keys())

    # --- Analysis helpers ---

    def _abundance_sql(self, level: str) -> str:
        """Build the long-form abundance SQL for a given level.

        Returns SQL that produces columns: sample_accession, feature_id, intensity.
        Dynamically detects whether the intensities struct uses 'label' (new
        schema) or 'channel' (old schema) so that old datasets keep working.
        """
        # The feature view still carries an intensities list<struct> that may use
        # 'label' (new) or 'channel' (old); detect for the peptide level. The pg
        # view is flattened since 1.1 (scalar p.label), so protein needs no probe.
        label_field = "label"
        if level == "peptide" and self.feature is not None:
            label_field = self.feature._intensity_label_field()

        if level == "protein":
            if self.pg is None or self.run is None:
                raise ValueError("level='protein' requires pg and run structures.")
            # Since QPX 1.1 the pg view is flattened: one row per label with scalar
            # p.label / p.intensity (no intensities list to UNNEST). Match each pg
            # row to the sample carrying its label. ID-only rows (null label /
            # intensity) never match a sample label and drop out naturally.
            return sql_build(
                """
            WITH numbered_pg AS MATERIALIZED (
                SELECT ROW_NUMBER() OVER () AS pg_row_id, *
                FROM pg
            ),
            pg_samples AS (
                SELECT DISTINCT p.pg_row_id,
                                rs.sample_accession,
                                rs.label AS intensity_label
                FROM numbered_pg p
                CROSS JOIN UNNEST(list_distinct(p.grouped_runs)) AS _g(run_file_name)
                JOIN run r USING (run_file_name)
                CROSS JOIN UNNEST(r.samples) AS _s(rs)
            )
            SELECT ps.sample_accession,
                   p.anchor_protein AS feature_id,
                   p.intensity
            FROM numbered_pg p
            JOIN pg_samples ps USING (pg_row_id)
            WHERE p.label = ps.intensity_label
              AND p.is_decoy = false
            """,
            )
        elif level == "peptide":
            if self.feature is None or self.run is None:
                raise ValueError("level='peptide' requires feature and run structures.")
            intensity_match = (
                "i.sample_accession = rs.sample_accession"
                if label_field == "channel"
                else "(i.label = rs.label OR len(r.samples) = 1)"
            )
            return sql_build(
                """
            SELECT rs.sample_accession,
                   f.sequence AS feature_id,
                   SUM(i.intensity) AS intensity
            FROM feature f,
                 run r,
                 UNNEST(r.samples) AS _t1(rs),
                 UNNEST(f.intensities) AS _t2(i)
            WHERE f.run_file_name = r.run_file_name
              AND $intensity_match
              AND f.is_decoy = false
            GROUP BY rs.sample_accession, f.sequence
            """,
                intensity_match=intensity_match,
            )
        else:
            raise ValueError(f"level must be 'protein' or 'peptide', got '{level}'")

    def intensity(self, level: str = "protein") -> QueryResult:
        """Lazy long-form intensity query — scalable for large datasets.

        Returns a QueryResult that stays lazy until you materialize it.
        DuckDB handles memory management, so this works on datasets that
        don't fit in memory.

        Parameters
        ----------
        level : str
            "protein" (uses pg + run) or "peptide" (uses feature + run).

        Returns
        -------
        QueryResult
            Long-form table with columns: sample_accession, feature_id, intensity.
            Call .to_df(), .to_arrow(), .to_polars(), or iterate row-by-row.

        Examples
        --------
        Lazy iteration (constant memory):
            for row in ds.intensity("protein"):
                sample, protein, intensity = row

        Materialize when data fits in memory:
            df = ds.intensity("protein").to_df()

        Write directly to Parquet (out-of-core):
            ds.intensity("protein").to_arrow()  # Arrow is more compact than pandas
        """
        sql = self._abundance_sql(level)
        return QueryResult(self._engine.execute(sql))

    def design_matrix(
        self,
        level: str = "protein",
        value_col: str = "intensity",
        fillna: float | None = 0.0,
        output_path: str | Path | None = None,
    ) -> "pd.DataFrame | Path":
        """Pivot intensity data into a samples-by-features matrix.

        For small-to-medium datasets, returns a pandas DataFrame.
        For large datasets, pass *output_path* to write the pivot directly
        to Parquet via DuckDB (out-of-core, no full pandas materialization).

        For streaming access to large data without pivoting, use
        :meth:`intensity` instead.

        Parameters
        ----------
        level : str
            "protein" (uses pg + run) or "peptide" (uses feature + run).
        value_col : str
            Which intensity column to pivot (default: "intensity").
        fillna : float | None
            Value to fill missing cells. Use *None* to keep NaN.
            Ignored when *output_path* is set.
        output_path : str | Path | None
            If set, writes the pivoted matrix directly to this Parquet path
            using DuckDB PIVOT (out-of-core). Returns the Path instead of a
            DataFrame.

        Returns
        -------
        pd.DataFrame | Path
            DataFrame (in-memory) or Path (when output_path is set).
        """
        import pandas as pd

        sql = self._abundance_sql(level)

        if output_path is not None:
            # Out-of-core: DuckDB PIVOT → Parquet, never touches pandas
            output_path = Path(output_path)
            pivot_sql = sql_build(
                """COPY (
                PIVOT ($base_sql) ON feature_id USING SUM(intensity)
            ) TO '$out_path' (FORMAT PARQUET)""",
                base_sql=sql,
                out_path=escape_path(str(output_path)),
            )
            self._engine.execute(pivot_sql)
            return output_path

        # In-memory path for small/medium datasets
        df = self._engine.execute(sql).fetchdf()

        if df.empty:
            return pd.DataFrame()

        matrix = df.pivot_table(
            index="sample_accession",
            columns="feature_id",
            values="intensity",
            aggfunc="sum",
        )
        matrix.columns.name = None

        if fillna is not None:
            matrix = matrix.fillna(fillna)

        return matrix

    # --- Validation ---
    def validate(self, structures: list[str] | None = None, *, strict: bool = False) -> dict[str, ValidationResult]:
        """Validate the dataset or specific structures against their schemas.

        Parameters
        ----------
        structures : list[str] | None
            Structure names to validate.  If *None*, validates all available.

        Returns
        -------
        dict[str, ValidationResult]
            Mapping of structure name to its validation result.
        """
        targets = structures or self.available_structures
        results: dict[str, ValidationResult] = {}
        for name in targets:
            struct = self._structures.get(name)
            if struct is None:
                result = ValidationResult(structure=name)
                result.issues.append(
                    ValidationIssue(
                        structure=name,
                        check="missing_structure",
                        severity="error",
                        column=None,
                        message=f"Structure '{name}' not found in dataset at {self.path}",
                    )
                )
                results[name] = result
            else:
                results[name] = struct.validate(strict=strict)

        # Cross-structure invariant: every pg.grouped_runs element must be a
        # real run.run_file_name. A grouped run that names no acquisition run is
        # dropped by every sample-joined view, so flag it as an error on pg.
        if "pg" in results and self.pg is not None and self.run is not None:
            self._check_grouped_runs_referential(results["pg"])
            self._check_grouped_runs_sample_mapping(results["pg"])

        return results

    def _check_grouped_runs_referential(self, pg_result: ValidationResult) -> None:
        """Flag pg.grouped_runs values that are not present in run.run_file_name.

        Appends an ``error`` issue to *pg_result* for each dangling grouped-run
        token. Any query failure (e.g. missing columns) is swallowed so schema
        validation is never masked by this referential check.
        """
        try:
            rows = self._engine.execute(
                """
                SELECT DISTINCT gr AS grouped_run
                FROM pg, UNNEST(pg.grouped_runs) AS _g(gr)
                WHERE gr IS NOT NULL
                  AND NOT EXISTS (SELECT 1 FROM run r WHERE r.run_file_name = gr)
                ORDER BY grouped_run
                """
            ).fetchall()
        except duckdb.Error as exc:
            # Do not fail silently: a query error here can mean genuine schema
            # drift (a referenced pg/run column was removed/renamed), which is how
            # the last release's bug hid. Surface it in the validation result.
            _log.warning("pg grouped_runs referential check could not run (possible schema drift): %s", exc)
            pg_result.issues.append(
                ValidationIssue(
                    structure="pg",
                    check="referential_check_skipped",
                    severity="warning",
                    column=None,
                    message=f"grouped_runs referential check could not run (possible schema drift): {exc}",
                )
            )
            return

        for (grouped_run,) in rows:
            pg_result.issues.append(
                ValidationIssue(
                    structure="pg",
                    check="dangling_grouped_run",
                    severity="error",
                    column="grouped_runs",
                    message=(
                        f"grouped_runs contains {grouped_run!r}, which is not a "
                        "run_file_name in run.parquet; PG rows keyed by it are "
                        "dropped from all sample-joined views"
                    ),
                )
            )

    def _check_grouped_runs_sample_mapping(self, pg_result: ValidationResult) -> None:
        """Require each PG intensity label to resolve to exactly one sample."""
        try:
            rows = self._engine.execute(
                """
                WITH numbered_pg AS MATERIALIZED (
                    SELECT ROW_NUMBER() OVER () AS pg_row_id, *
                    FROM pg
                ),
                pg_labels AS (
                    -- pg is flattened since 1.1: scalar p.label (no intensities
                    -- list). ID-only rows (null label) have no sample to resolve.
                    SELECT p.pg_row_id,
                           p.anchor_protein,
                           p.grouped_runs,
                           p.label
                    FROM numbered_pg p
                    WHERE p.label IS NOT NULL
                ),
                sample_matches AS (
                    SELECT pl.pg_row_id,
                           pl.label,
                           rs.sample_accession
                    FROM pg_labels pl
                    CROSS JOIN UNNEST(list_distinct(pl.grouped_runs)) AS _g(run_file_name)
                    JOIN run r USING (run_file_name)
                    CROSS JOIN UNNEST(r.samples) AS _s(rs)
                    WHERE rs.label = pl.label
                )
                SELECT pl.anchor_protein,
                       pl.grouped_runs,
                       pl.label,
                       COUNT(DISTINCT sm.sample_accession) AS sample_count,
                       LIST(DISTINCT sm.sample_accession)
                           FILTER (WHERE sm.sample_accession IS NOT NULL) AS samples
                FROM pg_labels pl
                LEFT JOIN sample_matches sm
                  ON pl.pg_row_id = sm.pg_row_id
                 AND pl.label = sm.label
                GROUP BY pl.pg_row_id,
                         pl.anchor_protein,
                         pl.grouped_runs,
                         pl.label
                HAVING COUNT(DISTINCT sm.sample_accession) != 1
                ORDER BY pl.anchor_protein, pl.label
                """
            ).fetchall()
        except duckdb.Error as exc:
            # Do not fail silently on drift (see _check_grouped_runs_referential).
            _log.warning("pg grouped_runs sample-mapping check could not run (possible schema drift): %s", exc)
            pg_result.issues.append(
                ValidationIssue(
                    structure="pg",
                    check="sample_mapping_check_skipped",
                    severity="warning",
                    column=None,
                    message=f"grouped_runs sample-mapping check could not run (possible schema drift): {exc}",
                )
            )
            return

        for anchor_protein, grouped_runs, label, sample_count, samples in rows:
            pg_result.issues.append(
                ValidationIssue(
                    structure="pg",
                    check="ambiguous_grouped_run_sample",
                    severity="error",
                    column="grouped_runs",
                    message=(
                        f"PG {anchor_protein!r} label {label!r} across "
                        f"grouped_runs={grouped_runs!r} resolves to "
                        f"{sample_count} samples ({samples or []!r}); exactly "
                        "one sample is required"
                    ),
                )
            )

    # --- Integrity ---
    def _require_local(self, operation: str) -> None:
        """Raise if the dataset is S3-backed; integrity and save require local paths."""
        if self._is_s3:
            raise NotImplementedError(
                f"{operation} is only supported for local datasets. This dataset is S3-backed (path={self.path!r})."
            )

    def compute_integrity(self) -> dict:
        """Compute checksums, row counts, and file sizes for all Parquet files.

        Only supported for local datasets. For S3-backed datasets, raises
        NotImplementedError.

        Returns
        -------
        dict
            Integrity fields suitable for writing to dataset.parquet (file_checksums,
            file_row_counts, file_sizes_bytes, total_structures, packaged_at).
        """
        self._require_local("compute_integrity")
        import hashlib
        from datetime import datetime, timezone

        import pyarrow.parquet as pq

        path = Path(self.path)
        checksums, row_counts, sizes = {}, {}, {}

        # Parquet files
        for f in sorted(path.glob("*.parquet")):
            name = f.name
            sizes[name] = f.stat().st_size
            sha = hashlib.sha256()
            with open(f, "rb") as fh:
                for chunk in iter(lambda: fh.read(65536), b""):
                    sha.update(chunk)
            checksums[name] = sha.hexdigest()
            try:
                row_counts[name] = pq.read_metadata(f).num_rows
            except Exception:
                row_counts[name] = -1

        # H5AD files (AnnData from downstream tools)
        for f in sorted(path.glob("*.h5ad")):
            name = f.name
            sizes[name] = f.stat().st_size
            sha = hashlib.sha256()
            with open(f, "rb") as fh:
                for chunk in iter(lambda: fh.read(65536), b""):
                    sha.update(chunk)
            checksums[name] = sha.hexdigest()
            try:
                import anndata

                row_counts[name] = anndata.read_h5ad(f, backed="r").n_obs
            except Exception:
                row_counts[name] = -1

        return {
            "file_checksums": checksums,
            "file_row_counts": row_counts,
            "file_sizes_bytes": sizes,
            "total_structures": len(checksums),
            "packaged_at": datetime.now(timezone.utc).isoformat(),
        }

    def verify_integrity(self) -> dict[str, list[str]]:
        """Verify dataset files against stored integrity data.

        Only supported for local datasets. For S3-backed datasets, returns
        a result with a single warning and no verification performed.

        Returns
        -------
        dict[str, list[str]]
            Dict with 'errors' and 'warnings' lists.
        """
        import hashlib

        errors, warnings = [], []
        if self._is_s3:
            warnings.append("Integrity verification is not supported for S3-backed datasets.")
            return {"errors": errors, "warnings": warnings}
        if self.dataset_meta is None:
            errors.append("No dataset.parquet found — cannot verify integrity")
            return {"errors": errors, "warnings": warnings}

        meta_df = self.dataset_meta.to_df()
        if meta_df.empty or "file_checksums" not in meta_df.columns:
            warnings.append("No integrity data stored in dataset.parquet")
            return {"errors": errors, "warnings": warnings}

        stored_checksums = meta_df["file_checksums"].iloc[0]
        if not isinstance(stored_checksums, dict):
            warnings.append("file_checksums is null")
            return {"errors": errors, "warnings": warnings}

        # Skip dataset.parquet itself — writing integrity changes its own checksum
        dataset_suffix = self._STRUCTURE_REGISTRY["dataset"][1]
        path = Path(self.path)
        for name, expected_sha in stored_checksums.items():
            if name.endswith(dataset_suffix):
                continue
            fpath = path / name
            if not fpath.exists():
                errors.append(f"Missing file: {name}")
                continue
            actual_sha = hashlib.sha256(fpath.read_bytes()).hexdigest()
            if actual_sha != expected_sha:
                errors.append(f"Checksum mismatch: {name}")

        return {"errors": errors, "warnings": warnings}

    # --- PRIDE enrichment ---
    def enrich_from_pride(self, project_accession: str | None = None) -> dict:
        """Fetch PRIDE metadata and update dataset.parquet.

        Retrieves project title, description, and PubMed ID from the
        PRIDE REST API and writes them into the dataset metadata.

        Parameters
        ----------
        project_accession : str or None
            PRIDE/ProteomeXchange accession. If *None*, reads from
            the existing dataset.parquet.

        Returns
        -------
        dict
            The fetched PRIDE metadata dict.

        Raises
        ------
        ValueError
            If no accession is provided and none is found in dataset.parquet,
            or if the accession is not found in PRIDE.
        """
        self._require_local("enrich_from_pride")

        # Resolve accession
        if project_accession is None:
            if self.dataset_meta is None:
                raise ValueError("No project_accession provided and no dataset.parquet found.")
            meta_df = self.dataset_meta.to_df()
            project_accession = meta_df["project_accession"].iloc[0]
            if not project_accession or project_accession == "unknown":
                raise ValueError("No valid project_accession in dataset.parquet. Provide one explicitly.")

        from qpx.core.pride import fetch_pride_metadata

        metadata = fetch_pride_metadata(project_accession)

        # Read existing dataset record, update fields, and rewrite
        if self.dataset_meta is not None:
            meta_df = self.dataset_meta.to_df()
            record = meta_df.iloc[0].to_dict()
        else:
            record = {
                "project_accession": project_accession,
                "creation_date": None,
            }

        record["project_title"] = metadata["project_title"]
        record["project_description"] = metadata["project_description"]
        record["pubmed_id"] = metadata["pubmed_id"]

        # Derive prefix from existing dataset file to overwrite in-place
        prefix = None
        if self.dataset_meta is not None:
            ds_suffix = self._STRUCTURE_REGISTRY["dataset"][1]
            existing_name = Path(self.dataset_meta._file_path).name
            if existing_name.endswith(ds_suffix):
                prefix = existing_name[: -len(ds_suffix)]

        self.save_structure([record], "dataset", prefix=prefix)
        self.refresh()

        return metadata

    # --- Write-back ---
    # Writer registry: name → (WriterClassName, file_suffix)
    _WRITER_REGISTRY = {
        "psm": ("PsmWriter", ".psm.parquet"),
        "feature": ("FeatureWriter", ".feature.parquet"),
        "pg": ("PgWriter", ".pg.parquet"),
        "mz": ("MzWriter", ".mz.parquet"),
        "sample": ("SampleWriter", ".sample.parquet"),
        "run": ("RunWriter", ".run.parquet"),
        "dataset": ("DatasetWriter", ".dataset.parquet"),
        "ontology": ("OntologyWriter", ".ontology.parquet"),
        "provenance": ("ProvenanceWriter", ".provenance.parquet"),
        "pepmap": ("PepMapWriter", ".pepmap.parquet"),
    }

    def save_structure(
        self,
        data,
        structure: str,
        prefix: str | None = None,
    ) -> Path:
        """Write validated Parquet back into the dataset directory.

        Only supported for local datasets. For S3-backed datasets, use the
        appropriate writer (e.g. ``FeatureWriter(path)``) with an explicit
        local output path.

        Parameters
        ----------
        data : list[dict] | pd.DataFrame | pa.Table
            Records to write. Schema-validated by the QPX writer.
        structure : str
            Structure name (e.g., "feature", "pg", "sample").
        prefix : str | None
            File name prefix. Defaults to the dataset directory name.

        Returns
        -------
        Path
            Path to the written Parquet file.
        """
        import pyarrow as pa

        self._require_local("save_structure")
        if structure not in self._WRITER_REGISTRY:
            raise ValueError(f"Unknown structure '{structure}'. Valid: {list(self._WRITER_REGISTRY.keys())}")

        writer_name, suffix = self._WRITER_REGISTRY[structure]

        import qpx.writers as writers_mod

        writer_cls = getattr(writers_mod, writer_name)

        path = Path(self.path)
        file_prefix = prefix or path.name
        output_path = path / f"{file_prefix}{suffix}"

        with writer_cls(output_path) as writer:
            if isinstance(data, list):
                writer.write_batch(data)
            elif isinstance(data, pa.Table):
                # Cast to the writer's schema to handle nullability/metadata differences
                casted = data.cast(writer.arrow_schema)
                writer.write_table(casted)
            else:
                writer.write_dataframe(data)

        return output_path

    # AnnData view type -> canonical file suffix (before .h5ad)
    _ANNDATA_VIEWS = {"ae", "de"}

    def save_anndata(
        self,
        adata,
        name: str | None = None,
        *,
        view: str | None = None,
    ) -> Path:
        """Write an AnnData object to the dataset directory as ``.h5ad``.

        Only supported for local datasets. For S3-backed datasets, write to
        a local path with ``adata.write_h5ad(path)``.

        Parameters
        ----------
        adata : anndata.AnnData
            The AnnData object to save.
        name : str, optional
            Base name for the file (without extension), e.g., "de_results".
            When omitted the dataset prefix is used.
        view : str, optional
            AnnData view type: ``"ae"`` (absolute expression) or ``"de"``
            (differential expression).  When provided the output file
            follows the QPX naming convention
            ``<prefix>.<view>.h5ad`` (e.g., ``PXD000000.ae.h5ad``).
            Ignored when *name* is explicitly given.

        Returns
        -------
        Path
            Path to the written .h5ad file.
        """
        self._require_local("save_anndata")
        path = Path(self.path)
        if name is not None:
            output_path = path / f"{name}.h5ad"
        elif view is not None:
            if view not in self._ANNDATA_VIEWS:
                raise ValueError(f"Unknown AnnData view '{view}'. Choose from: {sorted(self._ANNDATA_VIEWS)}")
            prefix = path.name
            output_path = path / f"{prefix}.{view}.h5ad"
        else:
            raise ValueError("Either 'name' or 'view' must be provided")
        adata.uns["qpx_version"] = QPX_SPEC_VERSION
        adata.uns["writer_version"] = __version__
        adata.write_h5ad(output_path)
        return output_path

    # --- External registration (for downstream tool integration) ---
    def register_external(self, name: str, file_path: str | Path) -> None:
        """Register an external Parquet file in the DuckDB engine."""
        self._engine.register_parquet(name, file_path)

    # --- Lifecycle ---
    def refresh(self) -> None:
        """Re-scan the dataset directory for new or updated files.

        Clears all cached view instances so next access gets fresh data.
        Call this after writing new structures or AnnData files.
        """
        for attr in [
            "_protein_view",
            "_peptide_view",
            "_identification_summary",
            "_run_summary",
            "_modification_view",
            "_qc_view",
            "_sample_summary",
            "_ae_view",
        ]:
            if hasattr(self, attr):
                delattr(self, attr)

        self._structures.clear()
        self._discover_and_register(None)

    def close(self):
        self._engine.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def __repr__(self):
        structs = ", ".join(self.available_structures)
        return f"Dataset('{self.path}', structures=[{structs}])"
