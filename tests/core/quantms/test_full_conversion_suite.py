#!/usr/bin/env python3
"""
Comprehensive Test Suite for QuantMS Conversions
Tests PSM, Feature, and Protein Group conversions for LFQ and TMT datasets
using the new unified CLI (qpx.cli.main).
"""

import os
import subprocess
import tempfile
import time
from pathlib import Path

import pyarrow.parquet as pq
import pytest

# New CLI entry point
CLI_MODULE = "qpx.cli.main"


def get_workspace_root():
    """Get the workspace root directory (tests/ for accessing examples/)"""
    return Path(__file__).parents[2]


def get_test_files():
    """Get paths to test files"""
    workspace_root = get_workspace_root()

    # Dataset paths
    lfq_dataset = workspace_root / "examples/quantms/dda-lfq-full"
    tmt_dataset = workspace_root / "examples/quantms/dda-plex-full"

    # File paths
    lfq_files = {
        "mztab": lfq_dataset / "PXD007683-LFQ.sdrf_openms_design_openms.mzTab.gz",
        "sdrf": lfq_dataset / "PXD007683-LFQ.sdrf.tsv",
        "msstats": lfq_dataset / "PXD007683-LFQ.sdrf_openms_design_msstats_in.csv.gz",
    }

    tmt_files = {
        "mztab": tmt_dataset / "PXD007683TMT.sdrf_openms_design_openms.mzTab.gz",
        "sdrf": tmt_dataset / "PXD007683-TMT.sdrf.tsv",
        "msstats": tmt_dataset / "PXD007683TMT.sdrf_openms_design_msstats_in.csv.gz",
    }

    return lfq_files, tmt_files


def run_command(cmd, description, workspace_root):
    """Run a command and return success/failure with timing"""
    print(f"\n[RUN] Running: {description}")
    print(f"Command: {' '.join(cmd)}")

    start_time = time.time()
    try:
        # Ensure UTF-8 encoding for subprocess
        env = os.environ.copy()
        env["PYTHONIOENCODING"] = "utf-8"

        result = subprocess.run(
            cmd,
            cwd=workspace_root,
            capture_output=True,
            text=True,
            env=env,
            encoding="utf-8",
            errors="replace",
            timeout=2400,
        )
        end_time = time.time()
        duration = end_time - start_time

        if result.returncode == 0:
            print(f"[OK] Success in {duration:.2f}s")
            return True, duration, result.stdout, result.stderr
        else:
            print(f"[FAIL] Failed in {duration:.2f}s")
            print(f"Error: {result.stderr}")
            return False, duration, result.stdout, result.stderr

    except subprocess.TimeoutExpired:
        print("[TIMEOUT] Timeout after 40 minutes")
        return False, 2400, "", "Timeout"
    except Exception as e:
        print(f"[ERROR] Exception: {str(e)}")
        return False, 0, "", str(e)


def analyze_output_file(file_path):
    """Analyze a parquet output file"""
    if not file_path.exists():
        return None

    try:
        df = pq.read_table(file_path).to_pandas()
        size_mb = file_path.stat().st_size / (1024 * 1024)

        return {
            "rows": len(df),
            "columns": len(df.columns),
            "size_mb": size_mb,
            "column_names": list(df.columns)[:10],  # First 10 columns
        }
    except Exception as e:
        return {"error": str(e)}


@pytest.mark.integration
def test_lfq_psm_conversion():
    """Test LFQ PSM conversion via unified quantms command"""
    workspace_root = get_workspace_root()
    lfq_files, _ = get_test_files()

    required_files = ["mztab", "sdrf"]
    for file_type in required_files:
        if not lfq_files[file_type].exists():
            pytest.skip(f"Test file not found: {lfq_files[file_type]}")

    with tempfile.TemporaryDirectory() as temp_dir:
        output_folder = Path(temp_dir)
        output_prefix = "LFQ_psm_test"

        cmd = [
            "python",
            "-m",
            CLI_MODULE,
            "convert",
            "quantms",
            "--mztab-path",
            str(lfq_files["mztab"]),
            "--sdrf-file",
            str(lfq_files["sdrf"]),
            "--output-folder",
            str(output_folder),
            "--output-prefix",
            output_prefix,
            "--structures",
            "psm",
            "--verbose",
        ]

        success, duration, stdout, stderr = run_command(cmd, "LFQ PSM Conversion", workspace_root)

        output_file = output_folder / f"{output_prefix}.psm.parquet"
        analysis = analyze_output_file(output_file)

        assert success, f"Command failed: {stderr}"
        assert output_file.exists(), f"Output file not found: {output_file}"
        assert analysis is not None, "Failed to analyze output file"
        assert "error" not in analysis, f"Analysis error: {analysis.get('error')}"
        assert analysis["rows"] > 0, "No data rows in output"


@pytest.mark.integration
def test_lfq_feature_conversion():
    """Test LFQ Feature conversion via unified quantms command"""
    workspace_root = get_workspace_root()
    lfq_files, _ = get_test_files()

    required_files = ["mztab", "sdrf", "msstats"]
    for file_type in required_files:
        if not lfq_files[file_type].exists():
            pytest.skip(f"Test file not found: {lfq_files[file_type]}")

    with tempfile.TemporaryDirectory() as temp_dir:
        output_folder = Path(temp_dir)
        output_prefix = "LFQ_feature_test"

        cmd = [
            "python",
            "-m",
            CLI_MODULE,
            "convert",
            "quantms",
            "--mztab-path",
            str(lfq_files["mztab"]),
            "--sdrf-file",
            str(lfq_files["sdrf"]),
            "--msstats-file",
            str(lfq_files["msstats"]),
            "--output-folder",
            str(output_folder),
            "--output-prefix",
            output_prefix,
            "--structures",
            "feature",
            "--verbose",
        ]

        success, duration, stdout, stderr = run_command(cmd, "LFQ Feature Conversion", workspace_root)

        output_file = output_folder / f"{output_prefix}.feature.parquet"
        analysis = analyze_output_file(output_file)

        assert success, f"Command failed: {stderr}"
        assert output_file.exists(), f"Output file not found: {output_file}"
        assert analysis is not None, "Failed to analyze output file"
        assert "error" not in analysis, f"Analysis error: {analysis.get('error')}"
        assert analysis["rows"] > 0, "No data rows in output"


@pytest.mark.integration
def test_lfq_protein_groups_conversion():
    """Test LFQ Protein Groups conversion via unified quantms command"""
    workspace_root = get_workspace_root()
    lfq_files, _ = get_test_files()

    required_files = ["mztab", "sdrf", "msstats"]
    for file_type in required_files:
        if not lfq_files[file_type].exists():
            pytest.skip(f"Test file not found: {lfq_files[file_type]}")

    with tempfile.TemporaryDirectory() as temp_dir:
        output_folder = Path(temp_dir)
        output_prefix = "LFQ_pg_test"

        cmd = [
            "python",
            "-m",
            CLI_MODULE,
            "convert",
            "quantms",
            "--mztab-path",
            str(lfq_files["mztab"]),
            "--sdrf-file",
            str(lfq_files["sdrf"]),
            "--msstats-file",
            str(lfq_files["msstats"]),
            "--output-folder",
            str(output_folder),
            "--output-prefix",
            output_prefix,
            "--structures",
            "feature,pg",
            "--verbose",
        ]

        success, duration, stdout, stderr = run_command(cmd, "LFQ Protein Groups Conversion", workspace_root)

        print(f"STDERR:\n{stderr}")

        output_file = output_folder / f"{output_prefix}.pg.parquet"
        analysis = analyze_output_file(output_file)

        assert success, f"Command failed: {stderr}"
        assert output_file.exists(), f"Output file not found: {output_file}"
        assert analysis is not None, "Failed to analyze output file"
        assert "error" not in analysis, f"Analysis error: {analysis.get('error')}"
        assert analysis["rows"] > 0, f"No data rows in output. Analysis: {analysis}"


@pytest.mark.integration
def test_tmt_psm_conversion():
    """Test TMT PSM conversion via unified quantms command"""
    workspace_root = get_workspace_root()
    _, tmt_files = get_test_files()

    required_files = ["mztab", "sdrf"]
    for file_type in required_files:
        if not tmt_files[file_type].exists():
            pytest.skip(f"Test file not found: {tmt_files[file_type]}")

    with tempfile.TemporaryDirectory() as temp_dir:
        output_folder = Path(temp_dir)
        output_prefix = "TMT_psm_test"

        cmd = [
            "python",
            "-m",
            CLI_MODULE,
            "convert",
            "quantms",
            "--mztab-path",
            str(tmt_files["mztab"]),
            "--sdrf-file",
            str(tmt_files["sdrf"]),
            "--output-folder",
            str(output_folder),
            "--output-prefix",
            output_prefix,
            "--structures",
            "psm",
            "--verbose",
        ]

        success, duration, stdout, stderr = run_command(cmd, "TMT PSM Conversion", workspace_root)

        output_file = output_folder / f"{output_prefix}.psm.parquet"
        analysis = analyze_output_file(output_file)

        assert success, f"Command failed: {stderr}"
        assert output_file.exists(), f"Output file not found: {output_file}"
        assert analysis is not None, "Failed to analyze output file"
        assert "error" not in analysis, f"Analysis error: {analysis.get('error')}"
        assert analysis["rows"] > 0, "No data rows in output"


@pytest.mark.integration
def test_tmt_feature_conversion():
    """Test TMT Feature conversion via unified quantms command"""
    workspace_root = get_workspace_root()
    _, tmt_files = get_test_files()

    required_files = ["mztab", "sdrf", "msstats"]
    for file_type in required_files:
        if not tmt_files[file_type].exists():
            pytest.skip(f"Test file not found: {tmt_files[file_type]}")

    with tempfile.TemporaryDirectory() as temp_dir:
        output_folder = Path(temp_dir)
        output_prefix = "TMT_feature_test"

        cmd = [
            "python",
            "-m",
            CLI_MODULE,
            "convert",
            "quantms",
            "--mztab-path",
            str(tmt_files["mztab"]),
            "--sdrf-file",
            str(tmt_files["sdrf"]),
            "--msstats-file",
            str(tmt_files["msstats"]),
            "--output-folder",
            str(output_folder),
            "--output-prefix",
            output_prefix,
            "--structures",
            "feature",
            "--verbose",
        ]

        success, duration, stdout, stderr = run_command(cmd, "TMT Feature Conversion", workspace_root)

        output_file = output_folder / f"{output_prefix}.feature.parquet"
        analysis = analyze_output_file(output_file)

        assert success, f"Command failed: {stderr}"
        assert output_file.exists(), f"Output file not found: {output_file}"
        assert analysis is not None, "Failed to analyze output file"
        assert "error" not in analysis, f"Analysis error: {analysis.get('error')}"
        assert analysis["rows"] > 0, "No data rows in output"


@pytest.mark.integration
def test_tmt_protein_groups_conversion():
    """Test TMT Protein Groups conversion via unified quantms command"""
    workspace_root = get_workspace_root()
    _, tmt_files = get_test_files()

    required_files = ["mztab", "sdrf", "msstats"]
    for file_type in required_files:
        if not tmt_files[file_type].exists():
            pytest.skip(f"Test file not found: {tmt_files[file_type]}")

    with tempfile.TemporaryDirectory() as temp_dir:
        output_folder = Path(temp_dir)
        output_prefix = "TMT_pg_test"

        cmd = [
            "python",
            "-m",
            CLI_MODULE,
            "convert",
            "quantms",
            "--mztab-path",
            str(tmt_files["mztab"]),
            "--sdrf-file",
            str(tmt_files["sdrf"]),
            "--msstats-file",
            str(tmt_files["msstats"]),
            "--output-folder",
            str(output_folder),
            "--output-prefix",
            output_prefix,
            "--structures",
            "feature,pg",
            "--verbose",
        ]

        success, duration, stdout, stderr = run_command(cmd, "TMT Protein Groups Conversion", workspace_root)

        print(f"STDERR:\n{stderr}")

        output_file = output_folder / f"{output_prefix}.pg.parquet"
        analysis = analyze_output_file(output_file)

        assert success, f"Command failed: {stderr}"
        assert output_file.exists(), f"Output file not found: {output_file}"
        assert analysis is not None, "Failed to analyze output file"
        assert "error" not in analysis, f"Analysis error: {analysis.get('error')}"
        assert analysis["rows"] > 0, f"No data rows in output. Analysis: {analysis}"


@pytest.mark.integration
def test_lfq_full_conversion():
    """Test LFQ full conversion (psm + feature + pg) in one command"""
    workspace_root = get_workspace_root()
    lfq_files, _ = get_test_files()

    required_files = ["mztab", "sdrf", "msstats"]
    for file_type in required_files:
        if not lfq_files[file_type].exists():
            pytest.skip(f"Test file not found: {lfq_files[file_type]}")

    with tempfile.TemporaryDirectory() as temp_dir:
        output_folder = Path(temp_dir)
        output_prefix = "LFQ_full_test"

        cmd = [
            "python",
            "-m",
            CLI_MODULE,
            "convert",
            "quantms",
            "--mztab-path",
            str(lfq_files["mztab"]),
            "--sdrf-file",
            str(lfq_files["sdrf"]),
            "--msstats-file",
            str(lfq_files["msstats"]),
            "--output-folder",
            str(output_folder),
            "--output-prefix",
            output_prefix,
            "--structures",
            "psm,feature,pg",
            "--verbose",
        ]

        success, duration, stdout, stderr = run_command(cmd, "LFQ Full Conversion", workspace_root)

        print(f"STDERR:\n{stderr}")

        assert success, f"Command failed: {stderr}"

        for view in ("psm", "feature", "pg"):
            output_file = output_folder / f"{output_prefix}.{view}.parquet"
            assert output_file.exists(), f"Missing output: {output_file}"
            analysis = analyze_output_file(output_file)
            assert analysis is not None, f"Failed to analyze {view} output"
            assert "error" not in analysis, f"{view} analysis error: {analysis.get('error')}"
            assert analysis["rows"] > 0, f"No data rows in {view} output"
