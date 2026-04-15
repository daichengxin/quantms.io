"""Shared base adapter for DIA-NN converters.

Provides common data-loading helpers used by both
:class:`DiannFeatureAdapter` and :class:`DiannPgAdapter`.
"""

from __future__ import annotations

from qpx.converters.base import BaseConverter
from qpx.core.sql import escape_path, sql_build


class DiaNNBaseAdapter(BaseConverter):
    """Base adapter with DIA-NN-specific loading utilities.

    Subclasses inherit ``_load_diann_report()`` so the logic is defined
    once rather than duplicated across the feature and PG adapters.
    """

    def _load_diann_report(self, path: str) -> None:
        """Create a DuckDB view over a DIA-NN report file.

        Uses ``CREATE VIEW`` so DuckDB reads from the file lazily with
        column pruning and predicate pushdown — avoiding loading the
        entire dataset into memory (critical for multi-GB reports).

        Supports both TSV (``report.tsv``) and Parquet (``report.parquet``)
        formats.  If the ``report`` view already exists it is skipped.

        Args:
            path: Filesystem path to the DIA-NN report file.
        """
        if self._table_exists("report"):
            self.logger.debug("report view already loaded -- skipping reload")
            return
        safe_path = escape_path(path)
        if path.endswith(".parquet"):
            self._conn.execute(
                sql_build(
                    "CREATE VIEW report AS SELECT * FROM read_parquet('$path')",
                    path=safe_path,
                )
            )
        else:
            self._conn.execute(
                sql_build(
                    """CREATE VIEW report AS
                SELECT * FROM read_csv_auto('$path',
                    delim='\t', header=true, auto_detect=true,
                    null_padding=true)""",
                    path=safe_path,
                )
            )
        count = self._conn.execute("SELECT COUNT(*) FROM report").fetchone()[0]
        self.logger.info(f"DIA-NN report view created ({count:,} rows)")
