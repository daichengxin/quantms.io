"""Tests for the pdc2qpx orchestration (download + CDAP + mz).

Downloads are skipped (``skip_download=True``); the CDAP ``.psm`` fixture is
staged into a ``<download_dir>/<study>/`` layout and, for the spectra case, a
mini mzML is generated with pyopenms. This exercises the orchestration wiring
without hitting the network.
"""

import shutil
from pathlib import Path

import pyarrow.parquet as pq
import pytest

_FIXTURE_CDAP = Path(__file__).resolve().parent.parent / "examples" / "cdap"


def _stage_study(tmp_path, study, with_mzml=False):
    """Build a ``download_dir/<study>/`` with fixture .psm (+ optional mini mzML)."""
    download_dir = tmp_path / "downloads"
    study_dir = download_dir / study
    study_dir.mkdir(parents=True)
    for psm in _FIXTURE_CDAP.glob("*.psm"):
        shutil.copy(psm, study_dir / psm.name)

    if with_mzml:
        pytest.importorskip("pyopenms")
        import pyopenms as oms

        exp = oms.MSExperiment()
        for ms_level, scan in [(1, 1), (2, 2)]:
            spectrum = oms.MSSpectrum()
            spectrum.setMSLevel(ms_level)
            spectrum.setRT(60.0 + scan)
            spectrum.setNativeID(f"controllerType=0 controllerNumber=1 scan={scan}")
            if ms_level >= 2:
                precursor = oms.Precursor()
                precursor.setMZ(500.0)
                precursor.setCharge(2)
                spectrum.setPrecursors([precursor])
            spectrum.set_peaks(([100.0, 200.0], [10.0, 20.0]))
            exp.addSpectrum(spectrum)
        oms.MzMLFile().store(str(study_dir / "run_x.mzML"), exp)

    return download_dir


def test_pdc2qpx_base(tmp_path):
    from qpx.pipeline.pdc2qpx import run_pdc2qpx

    study = "PDC_TEST"
    download_dir = _stage_study(tmp_path, study)
    out = tmp_path / "qpx"

    result = run_pdc2qpx(study, download_dir, out, skip_download=True)

    assert "mz" not in result
    feature = out / f"{study}.feature.parquet"
    assert feature.exists()
    assert pq.read_table(str(feature)).num_rows > 0


def test_pdc2qpx_with_spectra(tmp_path):
    from qpx.pipeline.pdc2qpx import run_pdc2qpx

    study = "PDC_TEST"
    download_dir = _stage_study(tmp_path, study, with_mzml=True)
    out = tmp_path / "qpx"

    result = run_pdc2qpx(study, download_dir, out, include_spectra=True, skip_download=True)

    assert result["mz"] == out / f"{study}.mz.parquet"
    mz = pq.read_table(str(result["mz"]))
    assert mz.num_rows == 2  # 1 MS1 + 1 MS2
    assert {"run_file_name", "scan"} <= set(mz.column_names)


def test_pdc2qpx_missing_psm_raises(tmp_path):
    from qpx.pipeline.pdc2qpx import run_pdc2qpx

    (tmp_path / "downloads" / "EMPTY").mkdir(parents=True)
    with pytest.raises(FileNotFoundError):
        run_pdc2qpx("EMPTY", tmp_path / "downloads", tmp_path / "qpx", skip_download=True)


def test_pdc2qpx_non_cdap_psm_raises(tmp_path):
    """A .psm whose header lacks CDAP signature columns is rejected early."""
    from qpx.pipeline.pdc2qpx import run_pdc2qpx

    study = "PDC_NONCDAP"
    study_dir = tmp_path / "downloads" / study
    study_dir.mkdir(parents=True)
    # A non-CDAP .psm: header without FileName/ScanNum/PeptideSequence/Protein.
    (study_dir / "bad.psm").write_text("colA\tcolB\tcolC\n1\t2\t3\n", encoding="utf-8")

    with pytest.raises(ValueError, match="not a CDAP-format"):
        run_pdc2qpx(study, tmp_path / "downloads", tmp_path / "qpx", skip_download=True)
