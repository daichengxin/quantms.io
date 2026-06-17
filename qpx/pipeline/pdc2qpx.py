"""pdc2qpx — one-shot PDC/CPTAC download + QPX conversion.

Given a PDC study accession, download its CDAP ``.psm`` files (and optionally the
mzML spectra), then convert them into a QPX dataset:

- ``.psm`` files  -> psm / feature / pg / dataset / ontology / provenance views
  (via :class:`qpx.converters.cdap.CdapConverter`)
- mzML files      -> a full-spectra ``mz.parquet``
  (via :meth:`qpx.transforms.spectra_mapping.SpectraMappingTransform.write_mz_parquet_from_dir`)

The download layer is :mod:`pridepy` (an optional dependency, ``qpx[pdc]``).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


def _download(study: str, file_type: str, download_dir: Path, download_threads: int) -> None:
    """Download one PDC file type into ``download_dir/<study>/`` via pridepy."""
    try:
        from pridepy.pdc.downloader import download_pdc_files
    except ImportError:
        raise ImportError("pridepy is required for pdc2qpx. Install it with: pip install qpx[pdc]") from None

    logger.info("Downloading PDC %s file_type=%s -> %s", study, file_type, download_dir)
    download_pdc_files(
        accession=study,
        file_type=file_type,
        output_folder=str(download_dir),
        skip_if_downloaded_already=True,
        checksum_check=True,
        download_threads=download_threads,
        retry=True,
    )


# Core CDAP columns that identify a CPTAC/PDC CDAP .psm header.
_CDAP_SIGNATURE_COLUMNS = ("FileName", "ScanNum", "PeptideSequence", "Protein")


def _assert_cdap_psm(study_dir: Path) -> None:
    """Pre-flight check that *study_dir* holds CDAP-format ``.psm`` files.

    Fails early with a clear message -- before the DuckDB conversion -- when the
    study provides no ``.psm`` files or when they are not CDAP output, instead of
    letting :class:`~qpx.converters.cdap.CdapConverter` fail later with an opaque
    "column not found" SQL error.
    """
    psm_files = sorted(study_dir.glob("*.psm")) if study_dir.is_dir() else []
    if not psm_files:
        raise FileNotFoundError(
            f"No .psm files found in {study_dir}. pdc2qpx needs CDAP .psm output; this PDC study may not provide harmonized PSMs."
        )
    with open(psm_files[0], encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    missing = [col for col in _CDAP_SIGNATURE_COLUMNS if col not in header]
    if missing:
        raise ValueError(
            f"{psm_files[0].name} is not a CDAP-format .psm (missing columns: "
            f"{missing}); pdc2qpx only supports CPTAC/PDC CDAP output."
        )


def run_pdc2qpx(
    study: str,
    download_dir: Path,
    output_folder: Path,
    *,
    include_spectra: bool = False,
    ms_levels: Optional[list[int]] = None,
    max_cpus: int = 24,
    max_memory: str = "16GB",
    download_threads: int = 24,
    skip_download: bool = False,
) -> dict[str, Path]:
    """Download a PDC study and convert it to a QPX dataset.

    Args:
        study: PDC study accession (e.g. ``"PDC000109"``).
        download_dir: Directory for downloaded files; PDC writes them under
            ``download_dir/<study>/``.
        output_folder: Directory for generated QPX parquet files.
        include_spectra: Also download mzML and produce a full-spectra
            ``<study>.mz.parquet``.
        ms_levels: MS levels for the mz view (``None`` = all). Precursor-level
            LFQ reanalysis needs MS1 as well as MS2.
        max_cpus: Threads for the CDAP conversion.
        max_memory: Memory limit for the CDAP conversion (e.g. ``"16GB"``).
        download_threads: Parallel HTTP Range threads per file for downloads.
        skip_download: Skip downloading and use files already present under
            ``download_dir/<study>/`` (useful for re-runs and tests).

    Returns:
        Mapping of produced outputs: ``{"base": output_folder}`` and, when
        *include_spectra* is set, ``{"mz": <study>.mz.parquet}``.
    """
    download_dir = Path(download_dir)
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    # PDC downloads land under <download_dir>/<study>/. The same directory holds
    # both .psm and .mzML.gz files; the CDAP and mz converters glob their own.
    study_dir = download_dir / study

    # 1. Base QPX from CDAP .psm files.
    if not skip_download:
        _download(study, "psm", download_dir, download_threads)
    _assert_cdap_psm(study_dir)

    from qpx.converters.cdap import CdapConverter

    logger.info("Converting CDAP .psm -> QPX base views (%s)", study)
    converter = CdapConverter(max_memory=max_memory, max_cpus=max_cpus)
    converter.convert(
        psm_dir=study_dir,
        output_folder=output_folder,
        output_prefix=study,
        project_accession=study,
    )
    outputs: dict[str, Path] = {"base": output_folder}

    # 2. Optional full-spectra mz.parquet from mzML files.
    if include_spectra:
        if not skip_download:
            _download(study, "mzml", download_dir, download_threads)
        if not any(p.name.lower().endswith((".mzml", ".mzml.gz")) for p in study_dir.iterdir()):
            raise FileNotFoundError(f"--include-spectra set but no mzML files found in {study_dir}")

        from qpx.transforms.spectra_mapping import SpectraMappingTransform

        mz_out = output_folder / f"{study}.mz.parquet"
        logger.info("Converting mzML -> %s", mz_out)
        SpectraMappingTransform(mzml_directory=study_dir).write_mz_parquet_from_dir(mz_out, ms_levels=ms_levels)
        outputs["mz"] = mz_out

    logger.info("pdc2qpx complete for %s -> %s", study, output_folder)
    return outputs
