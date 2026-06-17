"""``qpxc pdc2qpx`` — download a PDC/CPTAC study and convert it to QPX."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import click


@click.command("pdc2qpx")
@click.option("--study", required=True, help="PDC study accession (e.g. PDC000109)")
@click.option(
    "--download-dir",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="Directory for downloaded PDC files (written under <download-dir>/<study>/)",
)
@click.option(
    "--output-folder",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="Directory for generated QPX parquet files",
)
@click.option(
    "--include-spectra",
    is_flag=True,
    default=False,
    help="Also download mzML and produce a full-spectra <study>.mz.parquet",
)
@click.option(
    "--ms-levels",
    default=None,
    help="MS levels for the mz view (e.g. '2' or '1,2'). Default: all levels.",
)
@click.option("--max-cpus", default=24, show_default=True, help="Threads for the CDAP conversion")
@click.option("--max-memory", default="16GB", show_default=True, help="Memory limit for the CDAP conversion")
@click.option(
    "--download-threads",
    default=24,
    show_default=True,
    help="Parallel HTTP Range threads per file for downloads",
)
@click.option(
    "--skip-download",
    is_flag=True,
    default=False,
    help="Use files already present under <download-dir>/<study>/ (skip downloading)",
)
@click.option("--verbose", is_flag=True, help="Enable verbose logging")
def pdc2qpx_cmd(
    study: str,
    download_dir: Path,
    output_folder: Path,
    include_spectra: bool,
    ms_levels: Optional[str],
    max_cpus: int,
    max_memory: str,
    download_threads: int,
    skip_download: bool,
    verbose: bool,
):
    """Download a PDC/CPTAC study and convert it to a QPX dataset.

    Downloads the CDAP ``.psm`` files (and, with ``--include-spectra``, the mzML
    spectra) for one PDC study, then converts them: ``.psm`` to psm/feature/pg/...
    views and mzML to a full-spectra ``mz.parquet``. Requires the ``pdc`` extra
    (``pip install qpx[pdc]``) for the pridepy download layer.

    \b
    Examples:
        # Base QPX (psm/feature/pg) only
        qpxc pdc2qpx \\
            --study PDC000109 \\
            --download-dir ./downloads \\
            --output-folder ./qpx/PDC000109

        # Entire QPX including full spectra (for quantms reanalysis)
        qpxc pdc2qpx \\
            --study PDC000109 \\
            --download-dir ./downloads \\
            --output-folder ./qpx/PDC000109 \\
            --include-spectra --max-cpus 24
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    levels = [int(x.strip()) for x in ms_levels.split(",")] if ms_levels else None

    from qpx.pipeline.pdc2qpx import run_pdc2qpx

    outputs = run_pdc2qpx(
        study=study,
        download_dir=download_dir,
        output_folder=output_folder,
        include_spectra=include_spectra,
        ms_levels=levels,
        max_cpus=max_cpus,
        max_memory=max_memory,
        download_threads=download_threads,
        skip_download=skip_download,
    )

    click.echo(f"pdc2qpx complete for {study}. Output: {output_folder}")
    if "mz" in outputs:
        click.echo(f"Full spectra: {outputs['mz']}")
