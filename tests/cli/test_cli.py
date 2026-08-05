"""CLI command tests — convert, transform, query, info, validate, ontology."""

from click.testing import CliRunner

from qpx.cli.main import qpx_main

# ---------------------------------------------------------------------------
# Convert
# ---------------------------------------------------------------------------


def _assert_help(result, *options):
    """Verify CLI help renders and contains expected options."""
    if result.exit_code != 0:
        raise AssertionError(f"exit_code={result.exit_code}, output={result.output}")
    for opt in options:
        if opt not in result.output:
            raise AssertionError(f"Missing option {opt} in help output")


class TestDiaNNConvertCLI:
    def test_diann_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "diann", "--help"])
        _assert_help(result, "--report-path")


class TestMaxQuantConvertCLI:
    def test_maxquant_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "maxquant", "--help"])
        _assert_help(result, "--msms-file")


class TestFragPipeConvertCLI:
    def test_fragpipe_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "fragpipe", "--help"])
        _assert_help(result, "--psm-file")

    def test_fragpipe_help_shows_new_options(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "fragpipe", "--help"])
        _assert_help(
            result,
            "--ion-file",
            "--peptide-file",
            "--pg-file",
            "--experiment-annotation-file",
        )


class TestMzIdentMLConvertCLI:
    def test_mzidentml_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "mzidentml", "--help"])
        _assert_help(result, "--mzid-path")

    def test_mzidentml_help_shows_new_options(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "mzidentml", "--help"])
        _assert_help(result, "--mgf-path", "--include-spectra", "--project-accession")

    def test_mzidentml_help_shows_enrich_pride(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "mzidentml", "--help"])
        _assert_help(result, "--enrich-pride")


class TestSdrfConvertCLI:
    def test_sdrf_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "sdrf", "--help"])
        _assert_help(result, "--sdrf-file")


# ---------------------------------------------------------------------------
# Transform
# ---------------------------------------------------------------------------


class TestTransformGeneMapCLI:
    def test_genemap_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["transform", "gene-map", "--help"])
        _assert_help(result, "--fasta")


class TestTransformQuantifyCLI:
    def test_quantify_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["transform", "quantify", "--help"])
        _assert_help(result, "--feature-path", "--method", "--output")

    def test_quantify_help_shows_ibaq_options(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["transform", "quantify", "--help"])
        _assert_help(result, "--organism", "--ploidy", "--min-aa")


class TestTransformNormalizeAccessionsCLI:
    def test_normalize_accessions_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["transform", "normalize-accessions", "--help"])
        _assert_help(result, "--dataset", "--direction", "--fasta")

    def test_normalize_accessions_discovers_openms_prefix(self, tmp_path, monkeypatch):
        (tmp_path / "openms.feature.parquet").touch()
        (tmp_path / "openms.pg.parquet").touch()
        normalized = []

        def fake_normalize(**kwargs):
            normalized.append(kwargs["parquet_path"].name)
            return {"rows": 1, "accessions_changed": 0}

        monkeypatch.setattr("qpx.transforms.accession_normalizer.normalize_parquet", fake_normalize)

        runner = CliRunner()
        result = runner.invoke(
            qpx_main,
            ["transform", "normalize-accessions", "--dataset", str(tmp_path), "--in-place"],
        )

        assert result.exit_code == 0, result.output
        assert normalized == ["openms.feature.parquet", "openms.pg.parquet"]


class TestTransformUpdateMetadataCLI:
    def test_update_metadata_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["transform", "update-metadata", "--help"])
        _assert_help(result, "--dataset", "--sdrf", "--old-sdrf", "--force")


# ---------------------------------------------------------------------------
# Query
# ---------------------------------------------------------------------------


class TestQuerySqlCLI:
    def test_sql_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["query", "sql", "--help"])
        _assert_help(result, "--dataset-path", "--sql")


class TestQueryFilterCLI:
    def test_filter_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["query", "filter", "--help"])
        _assert_help(result, "--dataset-path", "--structure", "--condition")


class TestQueryHeadCLI:
    def test_head_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["query", "head", "--help"])
        _assert_help(result, "--dataset-path", "--structure")


# ---------------------------------------------------------------------------
# Info
# ---------------------------------------------------------------------------


class TestInfoCLI:
    def test_info_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["info", "--help"])
        _assert_help(result, "--dataset-path")

    def test_info_schema_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["info", "schema", "--help"])
        _assert_help(result, "--dataset-path")

    def test_info_metadata_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["info", "metadata", "--help"])
        _assert_help(result, "--file")


# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------


class TestValidateCLI:
    def test_validate_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["validate", "--help"])
        _assert_help(result, "--dataset-path", "--file", "--structure")


# ---------------------------------------------------------------------------
# Ontology
# ---------------------------------------------------------------------------


class TestOntologyCLI:
    def test_ontology_info_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["ontology", "info", "--help"])
        _assert_help(result, "--source")

    def test_ontology_search_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["ontology", "search", "--help"])
        _assert_help(result, "--source", "--top-k")

    def test_ontology_update_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["ontology", "update", "--help"])
        _assert_help(result, "--source")

    def test_ontology_build_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["ontology", "build", "--help"])
        _assert_help(result, "--source", "--all-sources")
