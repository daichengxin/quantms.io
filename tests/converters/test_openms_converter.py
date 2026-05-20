"""Tests for the OpenMS QPX converter (enrichment of -out_qpx output)."""

from __future__ import annotations

from pathlib import Path

import pyarrow.parquet as pq
import pytest

from qpx.converters.openms.converter import OpenMSConverter, _discover_parquet
from qpx.core.data.loader import load_schema
from qpx.writers.feature import FeatureWriter
from qpx.writers.pg import PgWriter
from qpx.writers.psm import PsmWriter
from tests.conftest import (
    make_feature_record,
    make_pg_record,
    make_psm_record,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_minimal_qpx(qpx_dir: Path, prefix: str = "quantms") -> None:
    """Write minimal QPX-compliant psm, feature, pg parquet files."""
    # PSM
    psm_records = [
        make_psm_record(sequence="PEPTIDEK", run_file_name="run_01", scan=[1001]),
        make_psm_record(sequence="ANOTHERK", run_file_name="run_01", scan=[1002], charge=3),
    ]
    with PsmWriter(qpx_dir / f"{prefix}.psm.parquet") as w:
        w.write_batch(psm_records)

    # Feature
    feat_records = [
        make_feature_record(sequence="PEPTIDEK", run_file_name="run_01", anchor_protein="P12345"),
        make_feature_record(sequence="ANOTHERK", run_file_name="run_01", anchor_protein="P12345", charge=3),
    ]
    with FeatureWriter(qpx_dir / f"{prefix}.feature.parquet") as w:
        w.write_batch(feat_records)

    # PG
    pg_records = [
        make_pg_record(anchor_protein="P12345", run_file_name="run_01"),
    ]
    with PgWriter(qpx_dir / f"{prefix}.pg.parquet") as w:
        w.write_batch(pg_records)


def _write_minimal_sdrf(sdrf_path: Path) -> None:
    """Write a minimal SDRF TSV for testing."""
    lines = [
        "source name\tcharacteristics[organism]\tcharacteristics[organism part]\t"
        "comment[data file]\tcomment[label]\tcomment[instrument]\t"
        "comment[cleavage agent details]\tcomment[fraction identifier]",
        "sample_1\tHomo sapiens\tliver\trun_01.raw\tAC=MS:1002038;NT=label free sample\t"
        "AC=MS:1000449;NT=LTQ Orbitrap\tAC=MS:1001251;NT=Trypsin\t1",
    ]
    sdrf_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Tests: _discover_parquet
# ---------------------------------------------------------------------------


class TestDiscoverParquet:
    def test_discovers_all_three(self, tmp_path):
        _write_minimal_qpx(tmp_path)
        found = _discover_parquet(tmp_path)
        assert set(found.keys()) == {"psm", "feature", "pg"}

    def test_empty_dir(self, tmp_path):
        found = _discover_parquet(tmp_path)
        assert not found

    def test_partial(self, tmp_path):
        """Only psm present."""
        psm_records = [make_psm_record()]
        with PsmWriter(tmp_path / "test.psm.parquet") as w:
            w.write_batch(psm_records)
        found = _discover_parquet(tmp_path)
        assert set(found.keys()) == {"psm"}


# ---------------------------------------------------------------------------
# Tests: OpenMSConverter
# ---------------------------------------------------------------------------


class TestOpenMSConverter:
    def test_convert_full_bundle(self, tmp_path):
        """Full enrichment: 3 core + SDRF -> 8 files."""
        qpx_dir = tmp_path / "openms_qpx"
        qpx_dir.mkdir()
        _write_minimal_qpx(qpx_dir)

        sdrf_path = tmp_path / "test.sdrf.tsv"
        _write_minimal_sdrf(sdrf_path)

        output = tmp_path / "output"
        converter = OpenMSConverter(qpx_dir=qpx_dir, sdrf_path=sdrf_path)
        converter.convert(output_folder=output, output_prefix="openms")

        # Check all 8 files exist
        expected_files = [
            "openms.psm.parquet",
            "openms.feature.parquet",
            "openms.pg.parquet",
            "openms.sample.parquet",
            "openms.run.parquet",
            "openms.ontology.parquet",
            "openms.provenance.parquet",
            "openms.dataset.parquet",
        ]
        for fname in expected_files:
            assert (output / fname).exists(), f"Missing: {fname}"

    def test_core_tables_validated(self, tmp_path):
        """Core tables should pass QPX schema validation after copy."""
        qpx_dir = tmp_path / "openms_qpx"
        qpx_dir.mkdir()
        _write_minimal_qpx(qpx_dir)

        sdrf_path = tmp_path / "test.sdrf.tsv"
        _write_minimal_sdrf(sdrf_path)

        output = tmp_path / "output"
        converter = OpenMSConverter(qpx_dir=qpx_dir, sdrf_path=sdrf_path)
        converter.convert(output_folder=output, output_prefix="openms")

        for view in ("psm", "feature", "pg"):
            schema = load_schema(view)
            table = pq.read_table(str(output / f"openms.{view}.parquet"))
            result = schema.validate_full(table)
            assert result.is_valid, f"{view} validation failed: {result.summary}"

    def test_provenance_structure(self, tmp_path):
        """Provenance should have OpenMS + QPX steps."""
        qpx_dir = tmp_path / "openms_qpx"
        qpx_dir.mkdir()
        _write_minimal_qpx(qpx_dir)

        sdrf_path = tmp_path / "test.sdrf.tsv"
        _write_minimal_sdrf(sdrf_path)

        output = tmp_path / "output"
        converter = OpenMSConverter(qpx_dir=qpx_dir, sdrf_path=sdrf_path)
        converter.convert(output_folder=output, output_prefix="openms")

        prov_table = pq.read_table(str(output / "openms.provenance.parquet"))
        assert prov_table.num_rows == 2
        tools = prov_table.column("tool_name").to_pylist()
        assert "OpenMS/ProteomicsLFQ" in tools
        assert "qpx" in tools

    def test_dataset_metadata(self, tmp_path):
        """Dataset should contain project accession."""
        qpx_dir = tmp_path / "openms_qpx"
        qpx_dir.mkdir()
        _write_minimal_qpx(qpx_dir)

        sdrf_path = tmp_path / "test.sdrf.tsv"
        _write_minimal_sdrf(sdrf_path)

        output = tmp_path / "output"
        converter = OpenMSConverter(qpx_dir=qpx_dir, sdrf_path=sdrf_path)
        converter.convert(
            output_folder=output,
            output_prefix="openms",
            project_accession="PXD001819",
        )

        ds_table = pq.read_table(str(output / "openms.dataset.parquet"))
        assert ds_table.num_rows == 1
        assert ds_table.column("project_accession").to_pylist()[0] == "PXD001819"
        assert ds_table.column("software_name").to_pylist()[0] == "OpenMS/ProteomicsLFQ"

    def test_no_qpx_files_raises(self, tmp_path):
        """Empty directory should raise FileNotFoundError."""
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        converter = OpenMSConverter(qpx_dir=empty_dir)
        with pytest.raises(FileNotFoundError, match="No QPX parquet files"):
            converter.convert(output_folder=tmp_path / "out")

    def test_same_dir_skips_copy(self, tmp_path):
        """When output_folder == qpx_dir, core files should not be copied."""
        _write_minimal_qpx(tmp_path, prefix="openms")

        sdrf_path = tmp_path / "test.sdrf.tsv"
        _write_minimal_sdrf(sdrf_path)

        converter = OpenMSConverter(qpx_dir=tmp_path, sdrf_path=sdrf_path)
        converter.convert(output_folder=tmp_path, output_prefix="openms")

        # Core files should still exist and metadata files should be added
        assert (tmp_path / "openms.psm.parquet").exists()
        assert (tmp_path / "openms.provenance.parquet").exists()

    def test_no_sdrf_skips_sample_run(self, tmp_path):
        """Without SDRF, sample/run should not be generated."""
        qpx_dir = tmp_path / "openms_qpx"
        qpx_dir.mkdir()
        _write_minimal_qpx(qpx_dir)

        output = tmp_path / "output"
        converter = OpenMSConverter(qpx_dir=qpx_dir)
        converter.convert(output_folder=output, output_prefix="openms")

        # Core + provenance + dataset + ontology should exist
        assert (output / "openms.psm.parquet").exists()
        assert (output / "openms.provenance.parquet").exists()
        assert (output / "openms.dataset.parquet").exists()
        # sample/run should NOT exist
        assert not (output / "openms.sample.parquet").exists()
        assert not (output / "openms.run.parquet").exists()

    def test_metadata_tables_validate(self, tmp_path):
        """All metadata tables should pass QPX schema validation."""
        qpx_dir = tmp_path / "openms_qpx"
        qpx_dir.mkdir()
        _write_minimal_qpx(qpx_dir)

        sdrf_path = tmp_path / "test.sdrf.tsv"
        _write_minimal_sdrf(sdrf_path)

        output = tmp_path / "output"
        converter = OpenMSConverter(qpx_dir=qpx_dir, sdrf_path=sdrf_path)
        converter.convert(output_folder=output, output_prefix="openms")

        for view in ("sample", "run", "ontology", "provenance", "dataset"):
            path = output / f"openms.{view}.parquet"
            if path.exists():
                schema = load_schema(view)
                table = pq.read_table(str(path))
                result = schema.validate_full(table)
                assert result.is_valid, f"{view} validation failed: {result.summary}"
