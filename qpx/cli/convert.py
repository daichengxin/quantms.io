"""
Convert subcommands — convert external tool outputs to QPX format.

Each command wraps a converter from ``qpx.converters.*``.

Subcommands:
    qpxc convert quantms    — QuantMS mzTab output to QPX
    qpxc convert diann      — DIA-NN report to QPX
    qpxc convert maxquant   — MaxQuant output to QPX
    qpxc convert fragpipe   — FragPipe output to QPX
    qpxc convert mzidentml  — mzIdentML (incl. XL-MS 1.3) to QPX
    qpxc convert sdrf       — SDRF to sample.parquet + run.parquet
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import click

logger = logging.getLogger("qpx.cli.convert")


def _maybe_enrich_pride(output_folder, project_accession: str | None, enrich: bool) -> None:
    """Optionally enrich a converted dataset with PRIDE metadata."""
    if not enrich:
        return
    if not project_accession:
        click.echo(
            "Warning: --enrich-pride requires --project-accession; skipping enrichment.",
            err=True,
        )
        return
    try:
        from qpx.dataset import Dataset

        ds = Dataset(output_folder)
        ds.enrich_from_pride(project_accession)
        ds.close()
        click.echo(f"Enriched dataset with PRIDE metadata for {project_accession}")
    except Exception as exc:
        click.echo(f"Warning: Could not enrich from PRIDE: {exc}", err=True)


# ---------------------------------------------------------------------------
# Convert group
# ---------------------------------------------------------------------------


@click.group()
def convert():
    """Convert external tool outputs to QPX format."""
    pass


# ---------------------------------------------------------------------------
# QuantMS
# ---------------------------------------------------------------------------


@convert.command("quantms")
@click.option(
    "--mztab-path",
    help="Input mzTab file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-file",
    help="SDRF metadata file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--msstats-file",
    help="MSstats input file path (required for feature and pg)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory for generated QPX files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--output-prefix",
    help="Prefix for output file names",
    default="quantms",
)
@click.option(
    "--structures",
    help="Comma-separated list of structures to produce (psm, feature, pg). Default: all.",
    default="psm,feature,pg",
)
@click.option(
    "--database-path",
    help="DuckDB database file path (reuse existing or create new)",
    type=click.Path(dir_okay=False, path_type=Path),
)
@click.option(
    "--project-accession",
    help="PRIDE / ProteomeXchange accession (e.g. PXD020192)",
)
@click.option(
    "--enrich-pride",
    help="Fetch project metadata from PRIDE API after conversion",
    is_flag=True,
    default=False,
)
@click.option(
    "--compression",
    type=click.Choice(["zstd", "snappy", "gzip", "none"], case_sensitive=False),
    default="zstd",
    show_default=True,
    help="Parquet compression codec.",
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_quantms_cmd(
    mztab_path: Path,
    sdrf_file: Path,
    msstats_file: Optional[Path],
    output_folder: Path,
    output_prefix: str,
    structures: str,
    database_path: Optional[Path],
    project_accession: Optional[str],
    enrich_pride: bool,
    compression: str,
    verbose: bool,
):
    """Convert QuantMS mzTab output to QPX format.

    Reads a QuantMS-produced mzTab file (with optional MSstats quantification)
    and writes QPX Parquet files for the requested data structures.

    \b
    Examples:
        # Convert everything (PSM + feature + protein groups)
        qpxc convert quantms \\
            --mztab-path data.mzTab \\
            --sdrf-file metadata.sdrf.tsv \\
            --msstats-file msstats_in.csv \\
            --output-folder ./qpx_output

        # Convert only PSMs
        qpxc convert quantms \\
            --mztab-path data.mzTab \\
            --sdrf-file metadata.sdrf.tsv \\
            --output-folder ./qpx_output \\
            --structures psm
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    from qpx.converters.quantms import QuantMSConverter

    converter = QuantMSConverter(
        mztab_path=mztab_path,
        sdrf_file=sdrf_file,
        msstats_file=msstats_file,
        database_path=database_path,
        compression=compression,
    )
    converter.convert(
        output_folder=output_folder,
        output_prefix=output_prefix,
        structures=[s.strip() for s in structures.split(",")],
        project_accession=project_accession,
    )

    _maybe_enrich_pride(output_folder, project_accession, enrich_pride)

    click.echo(f"QuantMS conversion complete. Output: {output_folder}")


# ---------------------------------------------------------------------------
# DIA-NN
# ---------------------------------------------------------------------------


@convert.command("diann")
@click.option(
    "--report-path",
    help="DIA-NN report file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-file",
    help="SDRF metadata file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--mzml-info-folder",
    help="Folder containing mzML info files (optional; scan/mz fields left empty if omitted)",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--qvalue-threshold",
    help="Q-value threshold for filtering",
    default=0.05,
    type=float,
)
@click.option(
    "--output-folder",
    help="Output directory for generated QPX files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option("--output-prefix", help="Prefix for output file names")
@click.option(
    "--pg-matrix-path",
    help="DIA-NN protein quantities matrix file (enables PG conversion)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--protein-file",
    help="Protein file for filtering",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--partitions",
    help="Field(s) for splitting output files (comma-separated)",
)
@click.option(
    "--duckdb-max-memory",
    help="Maximum memory for DuckDB engine (e.g., '4GB')",
)
@click.option(
    "--duckdb-threads",
    help="Number of threads for DuckDB engine",
    type=int,
)
@click.option(
    "--batch-size",
    help="Number of files to process simultaneously",
    default=100,
    type=int,
)
@click.option(
    "--standardized-intensities",
    help="Calculate standardized intensity metrics for PG output",
    is_flag=True,
    default=False,
)
@click.option(
    "--project-accession",
    help="PRIDE / ProteomeXchange accession (e.g. PXD020192)",
)
@click.option(
    "--enrich-pride",
    help="Fetch project metadata from PRIDE API after conversion",
    is_flag=True,
    default=False,
)
@click.option(
    "--compression",
    type=click.Choice(["zstd", "snappy", "gzip", "none"], case_sensitive=False),
    default="zstd",
    show_default=True,
    help="Parquet compression codec.",
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_diann_cmd(
    report_path: Path,
    sdrf_file: Path,
    mzml_info_folder: Path,
    qvalue_threshold: float,
    output_folder: Path,
    output_prefix: Optional[str],
    pg_matrix_path: Optional[Path],
    protein_file: Optional[Path],
    partitions: Optional[str],
    duckdb_max_memory: Optional[str],
    duckdb_threads: Optional[int],
    batch_size: int,
    standardized_intensities: bool,
    project_accession: Optional[str],
    enrich_pride: bool,
    compression: str,
    verbose: bool,
):
    """Convert DIA-NN report to QPX format.

    Reads a DIA-NN report.tsv file and converts feature-level quantification
    data into QPX Parquet format. When --pg-matrix-path is provided, also
    produces protein group output.

    \b
    Examples:
        # Feature conversion
        qpxc convert diann \\
            --report-path report.tsv \\
            --sdrf-file data.sdrf.tsv \\
            --mzml-info-folder ./mzml_info \\
            --output-folder ./qpx_output

        # Feature + protein groups
        qpxc convert diann \\
            --report-path report.tsv \\
            --sdrf-file data.sdrf.tsv \\
            --mzml-info-folder ./mzml_info \\
            --pg-matrix-path report.pg_matrix.tsv \\
            --output-folder ./qpx_output \\
            --standardized-intensities
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    output_folder.mkdir(parents=True, exist_ok=True)
    prefix = output_prefix or "diann"

    from qpx.converters.diann import DiaNNConverter

    converter = DiaNNConverter(
        report_path=report_path,
        sdrf_path=sdrf_file,
        duckdb_max_memory=duckdb_max_memory,
        duckdb_threads=duckdb_threads,
        compression=compression,
    )
    converter.convert_features(
        mzml_info_folder=mzml_info_folder,
        qvalue_threshold=qvalue_threshold,
        output_folder=output_folder,
        output_prefix=output_prefix,
        protein_file=protein_file,
        partitions=partitions,
        batch_size=batch_size,
    )
    if pg_matrix_path:
        converter.convert_pg(
            pg_matrix_path=pg_matrix_path,
            output_folder=output_folder,
            output_prefix=output_prefix,
            batch_size=batch_size,
            standardized_intensities=standardized_intensities,
        )

    converter.convert_sdrf(output_folder=output_folder, prefix=prefix)
    converter.write_ontology(output_folder, prefix=prefix)
    converter.write_provenance(output_folder, prefix=prefix)
    converter.write_dataset(output_folder, prefix=prefix, project_accession=project_accession)

    _maybe_enrich_pride(output_folder, project_accession, enrich_pride)

    click.echo(f"DIA-NN conversion complete. Output: {output_folder}")


# ---------------------------------------------------------------------------
# MaxQuant
# ---------------------------------------------------------------------------


@convert.command("maxquant")
@click.option(
    "--msms-file",
    help="MaxQuant msms.txt file (for PSM conversion)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--evidence-file",
    help="MaxQuant evidence.txt file (for feature conversion)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--protein-groups-file",
    help="MaxQuant proteinGroups.txt file (for PG conversion)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-file",
    help="SDRF metadata file (required for feature and PG)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory for generated QPX files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--output-prefix",
    help="Prefix for output file names",
)
@click.option(
    "--structures",
    help="Comma-separated list of structures to produce (psm, feature, pg). Default: all available.",
    default=None,
)
@click.option(
    "--protein-file",
    help="Protein list file for filtering feature output",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--batch-size",
    help="Processing batch size",
    default=100_000,
    type=int,
)
@click.option(
    "--n-workers",
    help="Number of parallel workers",
    type=int,
)
@click.option(
    "--memory-limit",
    help="Memory limit in GB",
    type=float,
)
@click.option(
    "--spectral-data",
    help="Include spectral data fields in PSM output",
    is_flag=True,
)
@click.option(
    "--standardized-intensities",
    help="Calculate standardized intensity metrics for PG output",
    is_flag=True,
    default=False,
)
@click.option(
    "--project-accession",
    help="PRIDE / ProteomeXchange accession (e.g. PXD020192)",
)
@click.option(
    "--enrich-pride",
    help="Fetch project metadata from PRIDE API after conversion",
    is_flag=True,
    default=False,
)
@click.option(
    "--compression",
    type=click.Choice(["zstd", "snappy", "gzip", "none"], case_sensitive=False),
    default="zstd",
    show_default=True,
    help="Parquet compression codec.",
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
@click.option(
    "--fixed-mod-only",
    help="Only include features with fixed modifications (Unmodified or Carbamidomethyl C)",
    is_flag=True,
    default=False,
)
def convert_maxquant_cmd(
    msms_file: Optional[Path],
    evidence_file: Optional[Path],
    protein_groups_file: Optional[Path],
    sdrf_file: Optional[Path],
    output_folder: Path,
    output_prefix: Optional[str],
    structures: Optional[str],
    protein_file: Optional[Path],
    batch_size: int,
    n_workers: Optional[int],
    memory_limit: Optional[float],
    spectral_data: bool,
    standardized_intensities: bool,
    project_accession: Optional[str],
    enrich_pride: bool,
    compression: str,
    verbose: bool,
    fixed_mod_only: bool,
):
    """Convert MaxQuant output to QPX format.

    Reads MaxQuant result files (msms.txt, evidence.txt, proteinGroups.txt)
    and writes corresponding QPX Parquet files.

    \b
    Examples:
        # Convert everything
        qpxc convert maxquant \\
            --msms-file msms.txt \\
            --evidence-file evidence.txt \\
            --protein-groups-file proteinGroups.txt \\
            --sdrf-file metadata.sdrf.tsv \\
            --output-folder ./qpx_output

        # Convert PSMs only
        qpxc convert maxquant \\
            --msms-file msms.txt \\
            --output-folder ./qpx_output \\
            --structures psm
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    output_folder.mkdir(parents=True, exist_ok=True)

    # Determine which structures to produce
    if structures:
        requested = {s.strip() for s in structures.split(",")}
    else:
        # Auto-detect from provided files
        requested = set()
        if msms_file:
            requested.add("psm")
        if evidence_file:
            requested.add("feature")
        if protein_groups_file:
            requested.add("pg")
    if not requested:
        raise click.ClickException(
            "No input files provided. Supply at least one of --msms-file, --evidence-file, or --protein-groups-file."
        )

    from qpx.converters.maxquant import MaxQuantConverter

    converter = MaxQuantConverter(memory_limit_gb=memory_limit, compression=compression)
    converter.convert(
        msms_file=msms_file,
        evidence_file=evidence_file,
        protein_groups_file=protein_groups_file,
        sdrf_file=sdrf_file,
        output_folder=output_folder,
        output_prefix=output_prefix,
        structures=list(requested),
        protein_file=protein_file,
        batch_size=batch_size,
        n_workers=n_workers,
        spectral_data=spectral_data,
        standardized_intensities=standardized_intensities,
        project_accession=project_accession,
        fixed_mod_only=fixed_mod_only,
    )

    _maybe_enrich_pride(output_folder, project_accession, enrich_pride)

    click.echo(f"MaxQuant conversion complete. Output: {output_folder}")


# ---------------------------------------------------------------------------
# FragPipe
# ---------------------------------------------------------------------------


@convert.command("fragpipe")
@click.option(
    "--psm-file",
    help="FragPipe psm.tsv file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--ion-file",
    help="FragPipe combined_ion.tsv file (for feature conversion)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--peptide-file",
    help="FragPipe combined_peptide.tsv file (for feature conversion)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--pg-file",
    help="FragPipe combined_protein.tsv file (for PG conversion)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-file",
    help="SDRF metadata file (for sample/run conversion)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory for generated QPX files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--output-prefix",
    help="Prefix for output file names",
)
@click.option(
    "--batch-size",
    help="Processing batch size",
    default=1_000_000,
    type=int,
)
@click.option(
    "--project-accession",
    help="PRIDE / ProteomeXchange accession (e.g. PXD020192)",
)
@click.option(
    "--enrich-pride",
    help="Fetch project metadata from PRIDE API after conversion",
    is_flag=True,
    default=False,
)
@click.option(
    "--compression",
    type=click.Choice(["zstd", "snappy", "gzip", "none"], case_sensitive=False),
    default="zstd",
    show_default=True,
    help="Parquet compression codec.",
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_fragpipe_cmd(
    psm_file: Optional[Path],
    ion_file: Optional[Path],
    peptide_file: Optional[Path],
    pg_file: Optional[Path],
    sdrf_file: Optional[Path],
    output_folder: Path,
    output_prefix: Optional[str],
    batch_size: int,
    project_accession: Optional[str],
    enrich_pride: bool,
    compression: str,
    verbose: bool,
):
    """Convert FragPipe output to QPX format.

    Reads FragPipe result files and converts them into QPX Parquet format.
    Supports psm.tsv, combined_ion.tsv, combined_peptide.tsv, and
    combined_protein.tsv.

    \b
    Examples:
        # Convert PSMs only
        qpxc convert fragpipe \\
            --psm-file psm.tsv \\
            --output-folder ./qpx_output

        # Convert features + protein groups
        qpxc convert fragpipe \\
            --ion-file combined_ion.tsv \\
            --pg-file combined_protein.tsv \\
            --sdrf-file metadata.sdrf.tsv \\
            --output-folder ./qpx_output
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if not any([psm_file, ion_file, peptide_file, pg_file]):
        raise click.ClickException(
            "No input files provided. Supply at least one of --psm-file, --ion-file, --peptide-file, or --pg-file."
        )

    output_folder.mkdir(parents=True, exist_ok=True)

    from qpx.converters.fragpipe import FragPipeConverter

    converter = FragPipeConverter(output_directory=output_folder, compression=compression)
    converter.convert(
        psm_file=psm_file,
        ion_file=ion_file,
        peptide_file=peptide_file,
        pg_file=pg_file,
        sdrf_file=sdrf_file,
        output_prefix=output_prefix,
        batch_size=batch_size,
        project_accession=project_accession,
    )

    _maybe_enrich_pride(output_folder, project_accession, enrich_pride)

    click.echo(f"FragPipe conversion complete. Output: {output_folder}")


# ---------------------------------------------------------------------------
# mzIdentML
# ---------------------------------------------------------------------------


@convert.command("mzidentml")
@click.option(
    "--mzid-path",
    help="Input mzIdentML (.mzid) file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory for generated QPX files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--output-prefix",
    help="Prefix for output file names",
    default="mzidentml",
)
@click.option(
    "--mgf-path",
    help="Optional MGF file for spectra attachment",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--include-spectra",
    help="Attach mz_array and intensity_array from MGF to PSM records",
    is_flag=True,
    default=False,
)
@click.option(
    "--project-accession",
    help="PRIDE / ProteomeXchange accession (e.g. PXD054720)",
)
@click.option(
    "--enrich-pride",
    help="Fetch project metadata from PRIDE API after conversion",
    is_flag=True,
    default=False,
)
@click.option(
    "--compression",
    type=click.Choice(["zstd", "snappy", "gzip", "none"], case_sensitive=False),
    default="zstd",
    show_default=True,
    help="Parquet compression codec.",
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_mzidentml_cmd(
    mzid_path: Path,
    output_folder: Path,
    output_prefix: str,
    mgf_path: Optional[Path],
    include_spectra: bool,
    project_accession: Optional[str],
    enrich_pride: bool,
    compression: str,
    verbose: bool,
):
    """Convert mzIdentML file to QPX dataset.

    Supports both standard mzIdentML (1.1/1.2) and mzIdentML 1.3 with
    cross-linking extensions (inter-peptide, looplinks, noncovalent).
    Produces a full QPX dataset including PSM, pepmap, provenance,
    ontology, and dataset metadata.

    \b
    Examples:
        # Convert a standard mzIdentML file
        qpxc convert mzidentml \\
            --mzid-path results.mzid \\
            --output-folder ./qpx_output

        # Convert an XL-MS mzIdentML 1.3 file
        qpxc convert mzidentml \\
            --mzid-path crosslinks.mzid \\
            --output-folder ./qpx_output \\
            --output-prefix xl_experiment

        # Convert with spectra from MGF
        qpxc convert mzidentml \\
            --mzid-path results.mzid \\
            --mgf-path spectra.mgf \\
            --include-spectra \\
            --output-folder ./qpx_output

        # Convert with project accession
        qpxc convert mzidentml \\
            --mzid-path results.mzid \\
            --output-folder ./qpx_output \\
            --project-accession PXD054720
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    from qpx.converters.mzidentml.converter import MzIdentMLConverter

    converter = MzIdentMLConverter(compression=compression)
    converter.convert(
        mzid_path=mzid_path,
        output_folder=output_folder,
        output_prefix=output_prefix,
        mgf_path=mgf_path,
        include_spectra=include_spectra,
        project_accession=project_accession,
    )

    _maybe_enrich_pride(output_folder, project_accession, enrich_pride)

    click.echo(f"mzIdentML conversion complete. Output: {output_folder}")


# ---------------------------------------------------------------------------
# SDRF
# ---------------------------------------------------------------------------


@convert.command("sdrf")
@click.option(
    "--sdrf-file",
    help="SDRF metadata file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory for generated QPX files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--output-prefix",
    help="Prefix for output file names",
    default="sdrf",
)
@click.option(
    "--compression",
    type=click.Choice(["zstd", "snappy", "gzip", "none"], case_sensitive=False),
    default="zstd",
    show_default=True,
    help="Parquet compression codec.",
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_sdrf_cmd(
    sdrf_file: Path,
    output_folder: Path,
    output_prefix: str,
    compression: str,
    verbose: bool,
):
    """Convert SDRF metadata to QPX sample.parquet and run.parquet.

    Reads a Sample and Data Relationship Format (SDRF) file and produces
    the QPX sample and run data structures as Parquet files.

    \b
    Example:
        qpxc convert sdrf \\
            --sdrf-file metadata.sdrf.tsv \\
            --output-folder ./qpx_output
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    output_folder.mkdir(parents=True, exist_ok=True)

    from qpx.converters.sdrf import SdrfConverter

    with SdrfConverter(compression=compression) as converter:
        converter.convert(
            sdrf_path=str(sdrf_file),
            sample_output=str(output_folder / f"{output_prefix}.sample.parquet"),
            run_output=str(output_folder / f"{output_prefix}.run.parquet"),
        )

    click.echo(f"SDRF conversion complete. Output: {output_folder}")
