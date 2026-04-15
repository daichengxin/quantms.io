"""Convert absolute-expression-2.0 quantms datasets into independent QPX datasets.

Each project folder on the PRIDE FTP contains quantms pipeline output
(mzTab, SDRF, MSstats CSV).  Unlike MSNet — where all projects merge into
one large QPX dataset — here every project becomes its own standalone QPX
dataset via the existing ``QuantMSConverter``.

Usage:
    # Convert all projects found in a local mirror
    python scripts/convert_absexpr.py --input /path/to/blood --output /path/to/output --all

    # Convert specific projects
    python scripts/convert_absexpr.py --input /path/to/blood --output /path/to/output \
        --projects PXD004682 PXD010154

    # Download from PRIDE FTP first, then convert
    python scripts/convert_absexpr.py --output /path/to/output --all --download
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
import time
from pathlib import Path

# Add QPX root to path so we can import qpx
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qpx.core.http import safe_urlopen

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
_log = logging.getLogger(__name__)

FTP_BASE = "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-collections/absolute-expression-2.0/blood"
DEFAULT_INPUT_DIR = "absexpr_blood"
DEFAULT_OUTPUT_DIR = "absexpr_qpx"


# ---------------------------------------------------------------------------
# File discovery helpers
# ---------------------------------------------------------------------------


def _find_single(directory: Path, glob_pattern: str) -> Path | None:
    """Return the single file matching *glob_pattern* in *directory*, or None."""
    matches = sorted(directory.glob(glob_pattern))
    if not matches:
        return None
    if len(matches) > 1:
        _log.warning("Multiple matches for %s in %s — using first", glob_pattern, directory)
    return matches[0]


def _discover_inputs(project_dir: Path) -> dict | None:
    """Locate the mzTab, SDRF, and MSstats files inside a project folder.

    Expected layout (quantms pipeline output):
        <project>/
            quant_tables/<accession>.mzTab
            quant_tables/<accession>_msstats_in.csv
            sdrf/<accession>.sdrf.tsv
    """
    quant_dir = project_dir / "quant_tables"
    sdrf_dir = project_dir / "sdrf"

    mztab = _find_single(quant_dir, "*.mzTab") if quant_dir.exists() else None
    if mztab is None:
        # Some older exports use mztab extensions
        mztab = _find_single(quant_dir, "*.mztab") if quant_dir.exists() else None

    sdrf = _find_single(sdrf_dir, "*.sdrf.tsv") if sdrf_dir.exists() else None

    if mztab is None or sdrf is None:
        _log.warning(
            "Skipping %s — missing %s",
            project_dir.name,
            "mzTab" if mztab is None else "SDRF",
        )
        return None

    msstats = _find_single(quant_dir, "*_msstats_in.csv") if quant_dir.exists() else None

    return {
        "mztab": mztab,
        "sdrf": sdrf,
        "msstats": msstats,
    }


def _list_projects(input_dir: Path) -> list[str]:
    """Return sorted list of project directory names."""
    return sorted(d.name for d in input_dir.iterdir() if d.is_dir() and not d.name.startswith("."))


# ---------------------------------------------------------------------------
# Download helper
# ---------------------------------------------------------------------------


def _download_file(url: str, dest: Path) -> bool:
    """Download a single file via HTTPS."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        _log.debug("Already exists: %s", dest.name)
        return True
    try:
        with safe_urlopen(url) as resp:
            dest.write_bytes(resp.read())
        return True
    except Exception:
        _log.debug("Failed to download %s", url, exc_info=True)
        return False


def _list_remote_files(url: str) -> list[str]:
    """Parse file names from an FTP directory index page."""
    import re

    try:
        with safe_urlopen(url) as resp:
            html = resp.read().decode("utf-8", errors="replace")
    except Exception:
        return []
    entries = re.findall(r'href="([^"/?][^"]*)"', html)
    return [e.rstrip("/") for e in entries if not e.startswith("?")]


def _download_project(project: str, input_dir: Path) -> bool:
    """Download the quant_tables and sdrf subdirs for a project using curl."""
    base_url = f"{FTP_BASE}/{project}"
    dest = input_dir / project
    _log.info("Downloading %s ...", project)

    ok = True
    for subdir in ("quant_tables", "sdrf"):
        dir_url = f"{base_url}/{subdir}/"
        files = _list_remote_files(dir_url)
        if not files:
            _log.warning("No files found in %s/%s", project, subdir)
            continue
        for fname in files:
            if fname.endswith("/"):
                continue  # skip sub-directories
            # Skip large files we don't need for conversion
            skip_ext = (".consensusXML", ".log", "_config.tsv", "_parsing.log")
            if any(fname.endswith(ext) for ext in skip_ext):
                _log.debug("Skipping %s (not needed)", fname)
                continue
            file_url = f"{dir_url}{fname}"
            file_dest = dest / subdir / fname
            if not _download_file(file_url, file_dest):
                _log.warning("Failed to download %s/%s/%s", project, subdir, fname)
                ok = False
            else:
                size_mb = file_dest.stat().st_size / (1024 * 1024)
                _log.info("  %s/%s (%.1f MB)", subdir, fname, size_mb)
    return ok


def _fetch_project_list() -> list[str]:
    """Fetch the list of project directories from the FTP index page."""
    import re

    url = f"{FTP_BASE}/"
    _log.info("Fetching project list from %s", url)
    try:
        with safe_urlopen(url) as resp:
            html = resp.read().decode("utf-8", errors="replace")
    except Exception as exc:
        _log.error("Could not fetch FTP index: %s", exc)
        return []

    # Parse directory links — PRIDE FTP index lists <a href="PXDxxxxxx/">
    dirs = re.findall(r'href="([A-Z0-9][^"/]+)/"', html)
    return sorted(set(dirs))


# ---------------------------------------------------------------------------
# Conversion
# ---------------------------------------------------------------------------


def _is_already_converted(output_dir: Path, project_name: str) -> bool:
    """Check if a project has already been converted (dataset.parquet exists)."""
    dataset_file = output_dir / project_name / "quantms.dataset.parquet"
    return dataset_file.exists()


def convert_project(
    project_name: str,
    inputs: dict,
    output_dir: Path,
) -> dict:
    """Convert a single project using QuantMSConverter.

    Returns a result dict with status and any error message.
    """
    from qpx.converters.quantms import QuantMSConverter

    project_output = output_dir / project_name
    t0 = time.perf_counter()
    result = {
        "project": project_name,
        "status": "success",
        "error": None,
        "mztab": str(inputs["mztab"]),
        "sdrf": str(inputs["sdrf"]),
        "msstats": str(inputs["msstats"]) if inputs["msstats"] else None,
        "output": str(project_output),
        "elapsed_s": 0.0,
    }

    structures = ["psm", "feature", "pg"] if inputs["msstats"] else ["psm"]

    try:
        converter = QuantMSConverter(
            mztab_path=inputs["mztab"],
            sdrf_file=inputs["sdrf"],
            msstats_file=inputs["msstats"],
        )
        converter.convert(
            output_folder=project_output,
            output_prefix="quantms",
            structures=structures,
            project_accession=project_name,
        )
        result["structures"] = structures
        result["elapsed_s"] = round(time.perf_counter() - t0, 1)
        _log.info("OK  %s → %s (%.1f s)", project_name, project_output, result["elapsed_s"])
    except Exception as exc:
        result["status"] = "error"
        result["error"] = str(exc)
        result["elapsed_s"] = round(time.perf_counter() - t0, 1)
        _log.error("FAIL %s — %s", project_name, exc)

    return result


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


def _write_report(results: list[dict], output_dir: Path) -> None:
    """Write a CSV summary report of all conversion results."""
    report_path = output_dir / "conversion_report.csv"
    fieldnames = ["project", "status", "elapsed_s", "structures", "mztab", "sdrf", "msstats", "output", "error"]

    with open(report_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for r in results:
            row = dict(r)
            if "structures" in row and isinstance(row["structures"], list):
                row["structures"] = ",".join(row["structures"])
            writer.writerow(row)

    ok = sum(1 for r in results if r["status"] == "success")
    fail = sum(1 for r in results if r["status"] == "error")
    skip = sum(1 for r in results if r["status"] == "skipped")
    _log.info("Report: %s  |  OK=%d  FAIL=%d  SKIP=%d", report_path, ok, fail, skip)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Convert absolute-expression-2.0/blood projects to QPX format.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=f"Root directory with project folders (default: {DEFAULT_INPUT_DIR})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output root directory (default: {DEFAULT_OUTPUT_DIR})",
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--all",
        action="store_true",
        dest="convert_all",
        help="Convert all projects found in the input directory",
    )
    group.add_argument(
        "--projects",
        nargs="+",
        help="List of project directory names to convert",
    )

    parser.add_argument(
        "--download",
        action="store_true",
        help="Download project data from PRIDE FTP before converting",
    )
    parser.add_argument(
        "--download-only",
        action="store_true",
        help="Only download data, do not convert",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        default=True,
        help="Skip projects that already have a dataset.parquet (default: True)",
    )
    parser.add_argument(
        "--no-skip-existing",
        action="store_false",
        dest="skip_existing",
        help="Re-convert projects even if output already exists",
    )
    parser.add_argument(
        "--structures",
        default="psm,feature,pg",
        help="Comma-separated list of structures to produce (default: psm,feature,pg)",
    )

    args = parser.parse_args()
    input_dir: Path = args.input.resolve()
    output_dir: Path = args.output.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine project list
    if args.download:
        projects = args.projects or _fetch_project_list()
        if not projects:
            _log.error("No projects to download")
            sys.exit(1)
        input_dir.mkdir(parents=True, exist_ok=True)
        for proj in projects:
            _download_project(proj, input_dir)
    elif args.convert_all:
        if not input_dir.exists():
            _log.error("Input directory does not exist: %s", input_dir)
            sys.exit(1)
        projects = _list_projects(input_dir)
    else:
        projects = args.projects

    if not projects:
        _log.error("No projects found in %s", input_dir)
        sys.exit(1)

    if args.download_only:
        _log.info("Download complete (--download-only). Exiting.")
        return

    _log.info("Projects to convert: %d", len(projects))

    results: list[dict] = []
    t0 = time.perf_counter()

    for i, proj in enumerate(projects, 1):
        _log.info("[%d/%d] Processing %s", i, len(projects), proj)
        proj_dir = input_dir / proj

        if not proj_dir.is_dir():
            _log.warning("Project directory not found: %s", proj_dir)
            results.append({"project": proj, "status": "skipped", "error": "directory not found"})
            continue

        if args.skip_existing and _is_already_converted(output_dir, proj):
            _log.info("Skipping %s — already converted (use --no-skip-existing to force)", proj)
            results.append({"project": proj, "status": "skipped", "error": "already converted"})
            continue

        inputs = _discover_inputs(proj_dir)
        if inputs is None:
            results.append({"project": proj, "status": "skipped", "error": "missing input files"})
            continue

        result = convert_project(proj, inputs, output_dir)
        results.append(result)

    elapsed = time.perf_counter() - t0
    _log.info("Finished %d projects in %.1f s", len(results), elapsed)
    _write_report(results, output_dir)


if __name__ == "__main__":
    main()
