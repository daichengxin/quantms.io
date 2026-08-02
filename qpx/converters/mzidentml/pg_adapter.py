"""mzIdentML PG adapter — converts ProteinDetectionList to pg.parquet.

Parses ``ProteinAmbiguityGroup`` elements from a mzIdentML file and
emits one ``pg.parquet`` row per protein group, using ``PgWriter``.
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path
from typing import Any

try:
    from lxml import etree
except ImportError:
    etree: Any = None

from qpx.converters.base import BaseConverter
from qpx.core.scores import normalize_score_name
from qpx.writers.pg import PgWriter

logger = logging.getLogger(__name__)

# CV accession constants
_CV_GROUP_REPRESENTATIVE = "MS:1002403"
_CV_LEADING_PROTEIN = "MS:1002401"
_CV_DISTINCT_PEPTIDES = "MS:1001097"
_CV_PROTEIN_FDR = "MS:1001364"

_DECOY_PREFIXES = ("rev_", "DECOY_", "decoy_")


class MzIdentMLPgAdapter(BaseConverter):
    """Convert mzIdentML ProteinDetectionList to QPX pg.parquet."""

    def convert(
        self,
        mzid_path: str | Path,
        output_path: str | Path,
        creator: str = "mzidentml",
    ) -> None:
        """Parse *mzid_path* and write QPX PG records to *output_path*.

        Args:
            mzid_path: Path to ``.mzid`` or ``.mzid.gz`` file.
            output_path: Destination ``pg.parquet`` path.
            creator: Creator tag written to Parquet metadata.
        """
        if etree is None:
            raise ImportError("lxml is required for mzIdentML conversion. Install it with: pip install qpx[mzidentml]")
        mzid_path = Path(mzid_path)
        run_file_name = mzid_path.name.replace(".mzid.gz", "").replace(".mzid", "")

        root, ns = self._parse_xml(mzid_path)
        dbseq_map = self._build_dbseq_map(root, ns)

        records = list(self._iter_protein_groups(root, ns, dbseq_map, run_file_name))

        if not records:
            logger.warning("No protein group records produced from %s", mzid_path)
            return

        with PgWriter(output_path, creator=creator, compression=self._compression) as writer:
            writer.write_batch(records)

        logger.info("Wrote %d protein groups to %s", len(records), output_path)

    # ------------------------------------------------------------------
    # XML parsing helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_xml(mzid_path: Path) -> tuple:
        """Parse mzIdentML XML, returns (root, namespace)."""
        if str(mzid_path).endswith(".gz"):
            with gzip.open(mzid_path, "rb") as fh:
                tree = etree.parse(fh)
        else:
            tree = etree.parse(str(mzid_path))
        root = tree.getroot()
        ns = root.nsmap.get(None, "")
        return root, ns

    @staticmethod
    def _build_dbseq_map(root, ns: str) -> dict[str, str]:
        """Build id → accession lookup from DBSequence elements."""
        return {d.get("id"): d.get("accession", d.get("id")) for d in root.iter(f"{{{ns}}}DBSequence")}

    def _iter_protein_groups(
        self,
        root,
        ns: str,
        dbseq_map: dict[str, str],
        run_file_name: str,
    ):
        """Yield one pg record dict per ProteinAmbiguityGroup."""
        for pag in root.iter(f"{{{ns}}}ProteinAmbiguityGroup"):
            record = self._build_pag_record(pag, ns, dbseq_map, run_file_name)
            if record is not None:
                yield record

    def _build_pag_record(self, pag, ns: str, dbseq_map: dict[str, str], run_file_name: str) -> dict | None:
        """Build a single pg record from a ProteinAmbiguityGroup element."""
        pdhs = pag.findall(f"{{{ns}}}ProteinDetectionHypothesis")
        if not pdhs:
            return None

        pg_accessions: list[str] = []
        anchor_protein: str | None = None
        peptides: list[dict] = []
        pg_qvalue: float | None = None
        additional_scores: list[dict] = []

        for pdh in pdhs:
            acc = dbseq_map.get(pdh.get("dBSequence_ref", ""), pdh.get("dBSequence_ref", ""))
            if acc:
                pg_accessions.append(acc)

            # Count PeptideHypothesis elements as peptide count for this protein
            pep_count = len(pdh.findall(f"{{{ns}}}PeptideHypothesis"))
            if acc:
                peptides.append({"protein_name": acc, "peptide_count": pep_count})

            # Extract cvParams
            for cv in pdh.findall(f"{{{ns}}}cvParam"):
                cv_acc = cv.get("accession", "")
                cv_name = cv.get("name", "")
                cv_val = cv.get("value")

                if cv_acc in (_CV_GROUP_REPRESENTATIVE, _CV_LEADING_PROTEIN):
                    if anchor_protein is None:
                        anchor_protein = acc
                elif cv_acc == _CV_PROTEIN_FDR:
                    try:
                        pg_qvalue = float(cv_val)
                    except (TypeError, ValueError):
                        pass
                elif cv_val is not None:
                    try:
                        score_val = float(cv_val)
                        additional_scores.append(
                            {
                                "score_name": normalize_score_name(cv_name or cv_acc),
                                "score_value": score_val,
                                "higher_better": None,
                            }
                        )
                    except (TypeError, ValueError):
                        pass

        if not pg_accessions:
            return None

        # Fallback anchor: first accession
        if anchor_protein is None:
            anchor_protein = pg_accessions[0]

        is_decoy = any(acc.startswith(_DECOY_PREFIXES) for acc in pg_accessions)

        return {
            "pg_accessions": pg_accessions,
            "pg_names": None,
            "gg_accessions": None,
            "gg_names": None,
            "gg_qvalue": None,
            "anchor_protein": anchor_protein,
            "grouped_runs": [run_file_name],
            "global_qvalue": None,
            "pg_qvalue": pg_qvalue,
            "intensities": None,
            "additional_intensities": None,
            "is_decoy": is_decoy,
            "contaminant": None,
            "peptides": peptides,
            "peptide_counts": None,
            "feature_counts": None,
            "sequence_coverage": None,
            "molecular_weight": None,
            "additional_scores": additional_scores if additional_scores else None,
            "cv_params": None,
        }
