"""Public ontology registry for QPX.

Provides bionty-style access to controlled vocabularies (PSI-MS, PRIDE CV,
and future ontologies) backed by Parquet files and DuckDB queries.

Usage::

    from qpx.core.ontology import PublicOntology

    ms = PublicOntology("psi_ms")
    ms.search("percolator")
    ms.lookup("MS:1001492")
    ms.validate(["percolator_score", "unknown"])
"""

from .registry import PublicOntology

__all__ = ["PublicOntology"]
