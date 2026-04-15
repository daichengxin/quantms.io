"""Spectra mapping transform — mzML spectral data extraction for PSM annotation.

This transform reads mzML files and extracts spectral arrays (m/z and intensity)
for identified scans, enabling PSM-level spectral annotation. It uses pyopenms
for mzML parsing with support for multiple native ID formats from different
instrument vendors.

The transform can:
1. Extract spectra for individual scans by scan number
2. Batch-extract spectra for all PSMs in a Dataset
3. Write annotated PSMs with spectral data to MZ Parquet files

Usage:
    mapper = SpectraMappingTransform(mzml_directory="/data/mzml/")

    # Extract a single spectrum
    spectrum = mapper.get_spectrum("run1", scan_number=12345)

    # Annotate PSMs in a Dataset with spectral arrays
    annotated_df = mapper.annotate_dataset_psms(dataset)

    # Write spectral data to .mz.parquet
    mapper.write_mz_parquet(dataset, "output.mz.parquet")
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# --- Native ID patterns for different instrument vendors ---
NATIVE_ID_PATTERNS = [
    "scan=(?<SCAN>\\d+)",  # Standard Thermo format
    "cycle=(?<SCAN>\\d+)",  # Waters/Agilent format
    "index=(?<SCAN>\\d+)",  # Generic index format
    "spectrum=(?<SCAN>\\d+)",  # Alternative format
]


class _MzMLCache:
    """
    Cache for loaded mzML experiments and spectrum lookup objects.

    Avoids reloading the same mzML file for multiple scan lookups.
    This is an internal implementation detail of SpectraMappingTransform.
    """

    def __init__(self):
        self._experiments: dict[str, object] = {}
        self._lookups: dict[str, object] = {}
        self._patterns: dict[str, str] = {}

    def get_spectrum(
        self,
        mzml_path: str,
        scan_number: int,
        native_id_pattern: Optional[str] = None,
    ) -> tuple[int, np.ndarray, np.ndarray]:
        """
        Get a spectrum from a cached mzML file.

        Args:
            mzml_path: Full path to the mzML file.
            scan_number: Scan number to retrieve.
            native_id_pattern: Optional custom regex for native ID parsing.

        Returns:
            Tuple of (num_peaks, mz_array, intensity_array).
            Returns (0, empty_array, empty_array) if scan not found.
        """
        try:
            import pyopenms as oms
            from pyopenms import SpectrumLookup
        except ImportError:
            raise ImportError("pyopenms is required for mzML parsing. Install it with: pip install pyopenms")

        # Load experiment if not cached
        if mzml_path not in self._experiments:
            logger.debug(f"Loading mzML file: {mzml_path}")
            exp = oms.MSExperiment()
            oms.MzMLFile().load(mzml_path, exp)
            self._experiments[mzml_path] = exp
            # Reset lookup for new file
            self._lookups.pop(mzml_path, None)
            self._patterns.pop(mzml_path, None)

        exp = self._experiments[mzml_path]

        # Build spectrum lookup if not cached
        if mzml_path not in self._lookups:
            patterns = [native_id_pattern] if native_id_pattern else NATIVE_ID_PATTERNS

            for pattern in patterns:
                try:
                    lookup = SpectrumLookup()
                    lookup.readSpectra(exp, pattern)
                    # Verify the pattern works by looking up the requested scan
                    lookup.findByScanNumber(scan_number)
                    self._lookups[mzml_path] = lookup
                    self._patterns[mzml_path] = pattern
                    logger.debug(f"Using native ID pattern: {pattern}")
                    break
                except Exception:
                    logger.debug("Native ID pattern %s failed for %s", pattern, mzml_path, exc_info=True)
                    continue

            if mzml_path not in self._lookups:
                logger.warning(f"Could not find working native ID pattern for {mzml_path}")
                return 0, np.array([], dtype=np.float32), np.array([], dtype=np.float32)

        lookup = self._lookups[mzml_path]

        try:
            index = lookup.findByScanNumber(scan_number)
        except (IndexError, RuntimeError):
            return 0, np.array([], dtype=np.float32), np.array([], dtype=np.float32)

        spectrum = exp.getSpectrum(index)
        mz_array, intensity_array = spectrum.get_peaks()

        return (
            len(mz_array),
            np.array(mz_array, dtype=np.float32),
            np.array(intensity_array, dtype=np.float32),
        )

    def clear(self):
        """Clear all cached experiments and lookups."""
        self._experiments.clear()
        self._lookups.clear()
        self._patterns.clear()


class SpectraMappingTransform:
    """
    Extract spectral arrays from mzML files and annotate PSM data.

    This transform maps scan numbers from PSM identifications to their
    corresponding m/z and intensity arrays in mzML files. It uses pyopenms
    for efficient mzML parsing with caching to avoid reloading files.

    The transform supports multiple instrument native ID formats (Thermo,
    Waters, Agilent, etc.) and automatically detects the correct format.

    Usage:
        mapper = SpectraMappingTransform(mzml_directory="/data/mzml/")

        # Get a single spectrum
        n_peaks, mz, intensity = mapper.get_spectrum("run1", 12345)

        # Annotate all PSMs in a dataset
        psm_df = mapper.annotate_dataset_psms(dataset)

        # Write MZ spectral data to Parquet
        mapper.write_mz_parquet(dataset, "output.mz.parquet")
    """

    # Common mzML file extensions to try
    MZML_EXTENSIONS = [".mzML", ".mzml", ".MZML"]

    def __init__(
        self,
        mzml_directory: Union[str, Path],
        native_id_pattern: Optional[str] = None,
    ):
        """
        Initialize the spectra mapping transform.

        Args:
            mzml_directory: Directory containing mzML files.
            native_id_pattern: Optional custom regex for native ID parsing.
                If not provided, common patterns are tried automatically.
        """
        self._mzml_dir = Path(mzml_directory)
        if not self._mzml_dir.is_dir():
            raise NotADirectoryError(f"mzML directory not found: {mzml_directory}")

        self._native_id_pattern = native_id_pattern
        self._cache = _MzMLCache()

    def _resolve_mzml_path(self, run_file_name: str) -> Optional[Path]:
        """
        Resolve a run file name to an mzML file path.

        Tries the run file name with common mzML extensions in the configured
        directory.

        Args:
            run_file_name: Run file name (without extension).

        Returns:
            Path to the mzML file, or None if not found.
        """
        for ext in self.MZML_EXTENSIONS:
            candidate = self._mzml_dir / f"{run_file_name}{ext}"
            if candidate.exists():
                return candidate

        # Also try the name as-is (might already have extension)
        candidate = self._mzml_dir / run_file_name
        if candidate.exists():
            return candidate

        return None

    def get_spectrum(
        self,
        run_file_name: str,
        scan_number: int,
    ) -> tuple[int, np.ndarray, np.ndarray]:
        """
        Get a single spectrum by run file name and scan number.

        Args:
            run_file_name: Name of the MS run (without mzML extension).
            scan_number: Scan number to retrieve.

        Returns:
            Tuple of (num_peaks, mz_array, intensity_array).
            Returns (0, empty_array, empty_array) if scan or file not found.
        """
        mzml_path = self._resolve_mzml_path(run_file_name)
        if mzml_path is None:
            logger.warning(f"mzML file not found for run: {run_file_name}")
            return 0, np.array([], dtype=np.float32), np.array([], dtype=np.float32)

        return self._cache.get_spectrum(
            str(mzml_path),
            scan_number,
            native_id_pattern=self._native_id_pattern,
        )

    def get_spectra_batch(
        self,
        run_file_name: str,
        scan_numbers: list[int],
    ) -> list[tuple[int, np.ndarray, np.ndarray]]:
        """
        Get multiple spectra from the same run file.

        More efficient than calling get_spectrum() repeatedly because the mzML
        file is loaded only once and kept in cache.

        Args:
            run_file_name: Name of the MS run.
            scan_numbers: List of scan numbers to retrieve.

        Returns:
            List of (num_peaks, mz_array, intensity_array) tuples.
        """
        mzml_path = self._resolve_mzml_path(run_file_name)
        if mzml_path is None:
            logger.warning(f"mzML file not found for run: {run_file_name}")
            empty = np.array([], dtype=np.float32)
            return [(0, empty.copy(), empty.copy()) for _ in scan_numbers]

        mzml_str = str(mzml_path)
        results = []
        for scan in scan_numbers:
            result = self._cache.get_spectrum(
                mzml_str,
                scan,
                native_id_pattern=self._native_id_pattern,
            )
            results.append(result)

        return results

    def annotate_dataset_psms(
        self,
        dataset,
        batch_size: int = 20,
    ) -> pd.DataFrame:
        """
        Annotate PSM data from a Dataset with spectral arrays.

        Retrieves m/z and intensity arrays for each PSM's scan number from
        the corresponding mzML file. Uses batched processing by run file
        to optimize mzML loading.

        Args:
            dataset: A qpx.Dataset with PSM data.
            batch_size: Number of run files to process in each batch.

        Returns:
            DataFrame with PSM data plus num_peaks, mz_array, intensity_array columns.
        """
        if dataset.psm is None:
            raise ValueError("Dataset does not contain PSM data.")

        annotated_chunks = []

        for run_files, batch_df in dataset.psm.iter_batches(
            partition_by="run_file_name",
            batch_size=batch_size,
        ):
            for run_file in run_files:
                run_mask = batch_df["run_file_name"] == run_file
                run_df = batch_df[run_mask].copy()

                if run_df.empty:
                    continue

                # Extract scan numbers from the scan column
                # In the new schema, scan is a list of int32 components
                scan_numbers = []
                for scan_val in run_df["scan"]:
                    if isinstance(scan_val, (list, np.ndarray)) and len(scan_val) > 0:
                        scan_numbers.append(int(scan_val[0]))
                    elif isinstance(scan_val, (int, np.integer)):
                        scan_numbers.append(int(scan_val))
                    else:
                        scan_numbers.append(0)

                # Batch-extract spectra for this run
                spectra = self.get_spectra_batch(run_file, scan_numbers)

                # Add spectral data to DataFrame
                run_df["num_peaks"] = [s[0] for s in spectra]
                run_df["mz_array"] = [s[1].tolist() for s in spectra]
                run_df["intensity_array"] = [s[2].tolist() for s in spectra]

                annotated_chunks.append(run_df)

        if not annotated_chunks:
            logger.warning("No PSMs were annotated with spectral data.")
            return pd.DataFrame()

        result = pd.concat(annotated_chunks, ignore_index=True)
        n_with_spectra = (result["num_peaks"] > 0).sum()
        logger.info(
            f"Annotated {n_with_spectra}/{len(result)} PSMs with spectral data ({n_with_spectra / len(result) * 100:.1f}%)"
        )
        return result

    def write_mz_parquet(
        self,
        dataset,
        output_path: Union[str, Path],
        ms_levels: Optional[list[int]] = None,
        batch_size: int = 10,
    ) -> Path:
        """
        Write spectral data from mzML files to .mz.parquet format.

        Reads mzML files for all runs in the Dataset and writes spectral data
        (m/z arrays, intensity arrays, precursor info, etc.) in the QPX MZ
        schema format.

        Args:
            dataset: A qpx.Dataset with run data (to identify mzML files).
            output_path: Path for the output .mz.parquet file.
            ms_levels: Optional list of MS levels to include (e.g., [1, 2]).
                If None, all MS levels are included.
            batch_size: Number of runs to process before flushing to disk.

        Returns:
            Path to the written file.
        """
        try:
            import pyopenms as oms
        except ImportError:
            raise ImportError("pyopenms is required for mzML parsing. Install it with: pip install pyopenms")

        from qpx.writers.mz import MzWriter

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        if dataset.run is None:
            raise ValueError("Dataset does not contain run data.")

        # Get run file names from the dataset
        run_df = dataset.run.to_df()
        run_file_names = run_df["run_file_name"].tolist()

        total_spectra = 0

        with MzWriter(output_path) as writer:
            for run_name in run_file_names:
                mzml_path = self._resolve_mzml_path(run_name)
                if mzml_path is None:
                    logger.warning(f"Skipping run {run_name}: mzML file not found")
                    continue

                logger.info(f"Processing mzML: {mzml_path}")
                exp = oms.MSExperiment()
                oms.MzMLFile().load(str(mzml_path), exp)

                records = []
                for i in range(exp.getNrSpectra()):
                    spectrum = exp.getSpectrum(i)
                    ms_level = spectrum.getMSLevel()

                    if ms_levels is not None and ms_level not in ms_levels:
                        continue

                    mz_array, intensity_array = spectrum.get_peaks()

                    # Build scan ID from native ID
                    native_id = spectrum.getNativeID()
                    scan_id = f"{run_name}:{native_id}" if native_id else f"{run_name}:index={i}"

                    # Extract precursor info for MS2+ scans
                    precursors = None
                    if ms_level >= 2 and spectrum.getPrecursors():
                        precursors = []
                        for precursor in spectrum.getPrecursors():
                            precursor_data = {
                                "selected_ion_mz": float(precursor.getMZ()),
                                "selected_ion_charge": (int(precursor.getCharge()) if precursor.getCharge() > 0 else None),
                                "selected_ion_intensity": (
                                    float(precursor.getIntensity()) if precursor.getIntensity() > 0 else None
                                ),
                                "isolation_window_target": float(precursor.getMZ()),
                                "isolation_window_lower": (
                                    float(precursor.getIsolationWindowLowerOffset())
                                    if precursor.getIsolationWindowLowerOffset() > 0
                                    else None
                                ),
                                "isolation_window_upper": (
                                    float(precursor.getIsolationWindowUpperOffset())
                                    if precursor.getIsolationWindowUpperOffset() > 0
                                    else None
                                ),
                                "spectrum_ref": None,
                            }
                            precursors.append(precursor_data)

                    record = {
                        "id": scan_id,
                        "ms_level": ms_level,
                        "centroid": spectrum.getType() == oms.SpectrumSettings.SpectrumType.CENTROID,
                        "scan_start_time": float(spectrum.getRT() / 60.0),  # Convert seconds to minutes
                        "inverse_ion_mobility": None,
                        "ion_injection_time": (
                            float(spectrum.getInstrumentSettings().getMetaValue("ion injection time"))
                            if spectrum.getInstrumentSettings().metaValueExists("ion injection time")
                            else 0.0
                        ),
                        "total_ion_current": (float(sum(intensity_array)) if len(intensity_array) > 0 else 0.0),
                        "precursors": precursors,
                        "mz": mz_array.tolist(),
                        "intensity": intensity_array.tolist(),
                        "cv_params": None,
                    }
                    records.append(record)

                    # Flush in batches
                    if len(records) >= 10000:
                        writer.write_batch(records)
                        total_spectra += len(records)
                        records = []

                # Flush remaining records for this run
                if records:
                    writer.write_batch(records)
                    total_spectra += len(records)

        logger.info(f"Wrote {total_spectra} spectra from {len(run_file_names)} runs to {output_path}")
        return output_path

    def close(self):
        """Release cached mzML data to free memory."""
        self._cache.clear()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
