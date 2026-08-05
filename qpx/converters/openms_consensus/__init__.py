"""Convert an OpenMS consensusXML (+ SDRF) directly to QPX.

Interim quantms production path while OpenMS `-out_qpx` is not yet emitting
format 1.1 (see PR bigbio/qpx#224). The consensusXML carries per-run peptide
feature intensities (`<element it=...>`), PSMs, and the protein-inference graph;
the SDRF supplies sample/label/fraction metadata and the ``grouped_runs`` units.
Protein-level abundance is intentionally NOT emitted here (the consensusXML has
no protein quantity; it lived only in the mzTab) — pg ships identification-only
until OpenMS `-out_qpx` provides the authoritative protein quant.
"""

from qpx.converters.openms_consensus.feature_adapter import consensus_features_to_records

__all__ = ["consensus_features_to_records"]
