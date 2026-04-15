"""Gene mapping transform — FASTA protein-to-gene annotation.

This transform maps protein accessions to gene names and gene accessions using
information parsed from UniProt FASTA headers and optionally from the MyGene.info
API for genomic accessions.

The FASTA header format expected (UniProt):
    >sp|P12345|PROT_HUMAN Some description GN=BRCA1 PE=1 SV=2

The 'GN=' field is extracted as the gene name (symbol).

Usage:
    mapping = GeneMappingTransform(fasta_path="proteins.fasta")

    # Annotate a Feature or PG structure
    annotated_df = mapping.annotate_dataframe(feature_df)

    # Annotate via Dataset (returns new DataFrame, does not modify originals)
    annotated_df = mapping.annotate_dataset_features(dataset)

    # Write annotated Parquet
    mapping.write_annotated_features(dataset, "output.feature.parquet")
"""

from __future__ import annotations

import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import Optional, Union

import pandas as pd

logger = logging.getLogger(__name__)


def _parse_gene_names_from_fasta(
    fasta_path: str,
    map_by: str = "accession",
) -> dict[str, set[str]]:
    """
    Parse FASTA file and build a protein -> gene name mapping.

    Extracts gene names from the 'GN=' field in UniProt FASTA headers.

    Args:
        fasta_path: Path to the FASTA file.
        map_by: How to key the mapping:
            - "accession": use UniProt accession (e.g., P12345)
            - "name": use UniProt entry name (e.g., PROT_HUMAN)

    Returns:
        Dict mapping protein identifier to set of gene names.
    """
    try:
        from Bio import SeqIO
    except ImportError:
        raise ImportError("Biopython is required for FASTA parsing. Install it with: pip install qpx[transforms]")

    gene_map: dict[str, set[str]] = defaultdict(set)

    for record in SeqIO.parse(fasta_path, "fasta"):
        # Parse the FASTA header: >db|ACCESSION|NAME description
        parts = record.id.split("|")

        if len(parts) >= 3:
            accession = parts[1]
            name = parts[2]
        elif len(parts) == 2:
            accession = parts[1]
            name = parts[1]
        else:
            accession = record.id
            name = record.id

        # Extract gene name from description using GN= field
        gene_list = re.findall(r"GN=(\S+)", record.description)
        gene_name = gene_list[0] if gene_list else None

        if map_by == "accession":
            gene_map[accession].add(gene_name)
        elif map_by == "name":
            gene_map[name].add(gene_name)
        else:
            # Map by both accession and name for maximum matching
            gene_map[accession].add(gene_name)
            gene_map[name].add(gene_name)

    logger.info(f"Parsed gene names for {len(gene_map)} protein identifiers from {fasta_path}")
    return gene_map


def _normalize_protein_list(accessions) -> Optional[list]:
    """Normalize a protein accessions value into a list suitable for _resolve_gene_names.

    - None / float NaN  → None
    - str               → [str]
    - dict / struct     → [item]   (list(dict) would give keys, not the struct)
    - iterable          → list(item)
    - other             → [item]
    """
    if accessions is None or (isinstance(accessions, float) and accessions != accessions):
        return None
    if isinstance(accessions, str):
        return [accessions]
    if isinstance(accessions, dict) or (
        hasattr(accessions, "dtype") and hasattr(accessions.dtype, "names") and accessions.dtype.names
    ):
        return [accessions]
    if hasattr(accessions, "__iter__"):
        return list(accessions)
    return [accessions]


def _extract_accession_from_raw(raw) -> Optional[str]:
    """Extract a bare accession string from a raw protein entry (str, dict, or struct).

    Returns None for None, NaN, or unsupported types.
    """
    if raw is None or (isinstance(raw, float) and raw != raw):
        return None
    if isinstance(raw, str):
        return raw
    if hasattr(raw, "__getitem__"):
        try:
            return raw["accession"]
        except (KeyError, TypeError):
            return None
    return None


def _resolve_gene_names(
    protein_accessions: list[str],
    gene_map: dict[str, set[str]],
) -> Optional[list[str]]:
    """
    Resolve gene names for a list of protein accessions.

    Args:
        protein_accessions: List of protein accessions from pg_accessions.
        gene_map: Protein-to-gene mapping from _parse_gene_names_from_fasta.

    Returns:
        List of gene names, or None if no mappings found.
    """
    if protein_accessions is None:
        return None

    gene_names = []
    for raw in protein_accessions:
        accession = _extract_accession_from_raw(raw)
        if accession is None:
            continue
        # Try full accession first, then extract short UniProt ID (sp|P12345|NAME → P12345)
        if accession not in gene_map:
            parts = accession.split("|")
            if len(parts) >= 2:
                accession = parts[1]
        if accession in gene_map:
            names = gene_map[accession]
            for name in names:
                if name is not None and name not in gene_names:
                    gene_names.append(name)

    return gene_names if gene_names else None


def _fetch_gene_accessions_from_api(
    gene_names: list[str],
    species: str = "human",
) -> dict[str, str]:
    """
    Fetch genomic accessions for gene symbols from MyGene.info API.

    Args:
        gene_names: List of gene symbols to look up.
        species: Species name for the query (default: "human").

    Returns:
        Dict mapping gene symbol to comma-separated genomic accession string.
    """
    if not gene_names:
        return {}

    try:
        import mygene
    except ImportError:
        logger.warning(
            "mygene package not installed. Gene accession lookup will be skipped. Install it with: pip install qpx[transforms]"
        )
        return {}

    mg = mygene.MyGeneInfo()

    try:
        results = mg.querymany(
            gene_names,
            scopes="symbol",
            species=species,
            fields="accession",
        )
    except Exception as e:
        logger.warning(f"MyGene.info API call failed: {e}")
        return {}

    accession_map: dict[str, str] = {}
    for result in results:
        query = result.get("query", "")
        accession_info = result.get("accession", {})
        if isinstance(accession_info, dict) and "genomic" in accession_info:
            genomic = accession_info["genomic"]
            if isinstance(genomic, list):
                accession_map[query] = ",".join(genomic)
            elif isinstance(genomic, str):
                accession_map[query] = genomic

    logger.info(f"Fetched genomic accessions for {len(accession_map)}/{len(gene_names)} genes")
    return accession_map


def _map_gene_accessions(
    gene_names: Optional[list[str]],
    accession_map: dict[str, str],
) -> Optional[list[str]]:
    """
    Map gene names to their genomic accessions.

    Args:
        gene_names: List of gene names.
        accession_map: Gene-to-accession mapping from _fetch_gene_accessions_from_api.

    Returns:
        List of genomic accessions, or None if no mappings found.
    """
    if gene_names is None:
        return None

    accessions = []
    for gene in gene_names:
        if gene in accession_map and accession_map[gene]:
            accessions.append(accession_map[gene])

    return accessions if accessions else None


class GeneMappingTransform:
    """
    Map protein accessions to gene names and gene accessions from FASTA.

    This transform reads UniProt FASTA headers to extract gene symbols (GN= field)
    and optionally queries MyGene.info for genomic accessions. It annotates
    QPX Feature or PG data structures with gg_names and gg_accessions columns.

    Usage:
        mapping = GeneMappingTransform(fasta_path="proteins.fasta")

        # Get the raw gene map
        gene_map = mapping.gene_map

        # Annotate a DataFrame (adds gg_names and gg_accessions columns)
        annotated_df = mapping.annotate_dataframe(df, protein_col="pg_accessions")

        # Annotate a Dataset's features
        feature_df = mapping.annotate_dataset_features(dataset)

        # Write annotated features to a new Parquet file
        mapping.write_annotated_features(dataset, "output.feature.parquet")
    """

    def __init__(
        self,
        fasta_path: Union[str, Path],
        map_by: str = "accession",
        species: str = "human",
        fetch_accessions: bool = True,
    ):
        """
        Initialize the gene mapping transform.

        Args:
            fasta_path: Path to the UniProt FASTA file.
            map_by: Mapping strategy ("accession" or "name").
            species: Species for MyGene.info lookup (default: "human").
            fetch_accessions: Whether to fetch genomic accessions from MyGene.info.
        """
        self._fasta_path = Path(fasta_path)
        if not self._fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

        self._map_by = map_by
        self._species = species
        self._fetch_accessions = fetch_accessions

        # Lazily computed
        self._gene_map: Optional[dict[str, set[str]]] = None
        self._accession_map: Optional[dict[str, str]] = None

    @property
    def gene_map(self) -> dict[str, set[str]]:
        """Protein-to-gene name mapping (lazy-loaded from FASTA)."""
        if self._gene_map is None:
            self._gene_map = _parse_gene_names_from_fasta(
                str(self._fasta_path),
                map_by=self._map_by,
            )
        return self._gene_map

    @property
    def accession_map(self) -> dict[str, str]:
        """Gene-to-genomic accession mapping (lazy-loaded from MyGene.info)."""
        if self._accession_map is None:
            if self._fetch_accessions:
                # Collect all unique gene names
                all_genes = set()
                for genes in self.gene_map.values():
                    for gene in genes:
                        if gene is not None:
                            all_genes.add(gene)
                self._accession_map = _fetch_gene_accessions_from_api(
                    sorted(all_genes),
                    species=self._species,
                )
            else:
                self._accession_map = {}
        return self._accession_map

    def annotate_dataframe(
        self,
        df: pd.DataFrame,
        protein_col: str = "pg_accessions",
        include_accessions: bool = True,
    ) -> pd.DataFrame:
        """
        Add gg_names and gg_accessions columns to a DataFrame.

        The protein_col should contain either:
        - A list of protein accessions (e.g., from QPX pg_accessions column)
        - A single protein accession string

        Args:
            df: Input DataFrame with protein identifiers.
            protein_col: Column name containing protein accessions.
            include_accessions: Whether to include gg_accessions from MyGene.info.

        Returns:
            DataFrame with added gg_names and gg_accessions columns.
        """
        if protein_col not in df.columns:
            raise ValueError(f"Column '{protein_col}' not found in DataFrame.")

        result = df.copy()
        gene_map = self.gene_map

        # Map protein accessions to gene names
        result["gg_names"] = result[protein_col].apply(
            lambda accessions: _resolve_gene_names(
                _normalize_protein_list(accessions),
                gene_map,
            )
        )

        # Map gene names to genomic accessions
        if include_accessions:
            acc_map = self.accession_map
            result["gg_accessions"] = result["gg_names"].apply(lambda names: _map_gene_accessions(names, acc_map))
        else:
            result["gg_accessions"] = None

        n_mapped = result["gg_names"].notna().sum()
        logger.info(f"Mapped gene names for {n_mapped}/{len(result)} rows ({n_mapped / len(result) * 100:.1f}%)")
        return result

    def annotate_dataset_features(
        self,
        dataset,
        include_accessions: bool = True,
    ) -> pd.DataFrame:
        """
        Annotate a Dataset's Feature data with gene names and accessions.

        Materializes the Feature data as a DataFrame, adds gene annotations,
        and returns the annotated DataFrame.

        Args:
            dataset: A qpx.Dataset with feature data.
            include_accessions: Whether to include gg_accessions from MyGene.info.

        Returns:
            Annotated Feature DataFrame with gg_names and gg_accessions.
        """
        if dataset.feature is None:
            raise ValueError("Dataset does not contain feature data.")

        feature_df = dataset.feature.to_df()
        return self.annotate_dataframe(
            feature_df,
            protein_col="pg_accessions",
            include_accessions=include_accessions,
        )

    def annotate_dataset_pg(
        self,
        dataset,
        include_accessions: bool = True,
    ) -> pd.DataFrame:
        """
        Annotate a Dataset's PG data with gene names and accessions.

        Args:
            dataset: A qpx.Dataset with PG data.
            include_accessions: Whether to include gg_accessions from MyGene.info.

        Returns:
            Annotated PG DataFrame with gg_names and gg_accessions.
        """
        if dataset.pg is None:
            raise ValueError("Dataset does not contain PG data.")

        pg_df = dataset.pg.to_df()
        return self.annotate_dataframe(
            pg_df,
            protein_col="pg_accessions",
            include_accessions=include_accessions,
        )

    def write_annotated_features(
        self,
        dataset,
        output_path: Union[str, Path],
        include_accessions: bool = True,
    ) -> Path:
        """
        Write gene-annotated Feature data to a new Parquet file.

        Uses the FeatureWriter to produce a schema-validated output file.

        Args:
            dataset: A qpx.Dataset with feature data.
            output_path: Path for the output .feature.parquet file.
            include_accessions: Whether to include gg_accessions from MyGene.info.

        Returns:
            Path to the written file.
        """
        from qpx.writers.feature import FeatureWriter

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        annotated_df = self.annotate_dataset_features(
            dataset,
            include_accessions=include_accessions,
        )

        with FeatureWriter(output_path) as writer:
            writer.write_dataframe(annotated_df)

        logger.info(f"Wrote gene-annotated features to {output_path}")
        return output_path

    def write_annotated_pg(
        self,
        dataset,
        output_path: Union[str, Path],
        include_accessions: bool = True,
    ) -> Path:
        """
        Write gene-annotated PG data to a new Parquet file.

        Uses the PgWriter to produce a schema-validated output file.

        Args:
            dataset: A qpx.Dataset with PG data.
            output_path: Path for the output .pg.parquet file.
            include_accessions: Whether to include gg_accessions from MyGene.info.

        Returns:
            Path to the written file.
        """
        from qpx.writers.pg import PgWriter

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        annotated_df = self.annotate_dataset_pg(
            dataset,
            include_accessions=include_accessions,
        )

        with PgWriter(output_path) as writer:
            writer.write_dataframe(annotated_df)

        logger.info(f"Wrote gene-annotated PG to {output_path}")
        return output_path

    def get_protein_sequences(self) -> dict[str, str]:
        """
        Parse protein sequences from FASTA for position mapping.

        Returns:
            Dict mapping protein accession to amino acid sequence.
        """
        try:
            from Bio import SeqIO
        except ImportError:
            raise ImportError("Biopython is required for FASTA parsing. Install it with: pip install biopython")

        protein_dict: dict[str, str] = {}
        for record in SeqIO.parse(str(self._fasta_path), "fasta"):
            parts = record.id.split("|")
            accession = parts[1] if len(parts) >= 2 else record.id
            protein_dict[accession] = str(record.seq)

        return protein_dict
