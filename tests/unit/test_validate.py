"""Tests for the QPX validation feature.

Covers: ValidationIssue, ValidationResult dataclasses,
ViewSchema.validate_full(), BaseStructure.validate(),
Dataset.validate(), and the CLI validate command.
"""

import pyarrow as pa
from click.testing import CliRunner

from qpx.core.data import (
    DatasetSchema,
    FeatureSchema,
    ValidationIssue,
    ValidationResult,
)


def test_validation_result_dataclass():
    """ValidationResult: is_valid, warnings-only, errors, errors/warnings properties."""
    # No issues = valid
    r = ValidationResult(structure="feature")
    assert r.is_valid is True

    # Warnings only = still valid
    r = ValidationResult(
        structure="feature",
        issues=[ValidationIssue("feature", "null_values", "warning", "col", "msg")],
    )
    assert r.is_valid is True

    # Error = invalid
    r = ValidationResult(
        structure="feature",
        issues=[ValidationIssue("feature", "missing_column", "error", "col", "msg")],
    )
    assert r.is_valid is False

    # Errors/warnings filter
    r = ValidationResult(
        structure="feature",
        issues=[
            ValidationIssue("feature", "missing_column", "error", "a", "err"),
            ValidationIssue("feature", "null_values", "warning", "b", "warn"),
            ValidationIssue("feature", "type_mismatch", "error", "c", "err2"),
        ],
    )
    assert len(r.errors) == 2
    assert len(r.warnings) == 1


def test_validate_full():
    """validate_full: valid, missing column, type mismatch, optional absent, nulls, duplicate PK."""
    schema = FeatureSchema.get_arrow_schema()

    # Valid table
    arrays = {f.name: pa.nulls(1, type=f.type) for f in schema}
    table = pa.table(arrays, schema=schema)
    result = FeatureSchema.validate_full(table)
    assert result.is_valid

    # Missing required column
    fields_to_keep = [f for f in schema if f.name != "sequence"]
    arrays = {f.name: pa.nulls(1, type=f.type) for f in fields_to_keep}
    table = pa.table(arrays, schema=pa.schema(fields_to_keep))
    result = FeatureSchema.validate_full(table)
    assert not result.is_valid
    assert any(i.check == "missing_column" and "sequence" in i.message for i in result.issues)

    # Type mismatch
    wrong_fields = [pa.field(f.name, pa.string()) if f.name == "charge" else f for f in schema]
    arrays = {f.name: pa.nulls(1, type=f.type) for f in wrong_fields}
    arrays["charge"] = pa.array(["wrong"], type=pa.string())
    table = pa.table(arrays, schema=pa.schema(wrong_fields))
    result = FeatureSchema.validate_full(table)
    assert not result.is_valid
    assert any(i.check == "type_mismatch" and "charge" in i.message for i in result.issues)

    # Optional column absent
    fields_to_keep = [f for f in schema if f.name != "pg_global_qvalue"]
    arrays = {f.name: pa.nulls(1, type=f.type) for f in fields_to_keep}
    table = pa.table(arrays, schema=pa.schema(fields_to_keep))
    result = FeatureSchema.validate_full(table)
    assert not any(i.check == "missing_column" and "pg_global_qvalue" in i.message for i in result.issues)

    # Null in non-nullable column = warning
    arrays = {f.name: pa.nulls(1, type=f.type) for f in schema}
    arrays["sequence"] = pa.array([None], type=pa.string())
    table = pa.table(arrays, schema=schema)
    result = FeatureSchema.validate_full(table)
    null_issues = [i for i in result.issues if i.check == "null_values" and i.column == "sequence"]
    assert len(null_issues) == 1
    assert null_issues[0].severity == "warning"

    # Duplicate PK = warning
    ds_schema = DatasetSchema.get_arrow_schema()
    arrays = {}
    for f in ds_schema:
        if f.name == "project_accession":
            arrays[f.name] = pa.array(["PXD001", "PXD001"], type=f.type)
        elif f.name == "creation_date":
            arrays[f.name] = pa.array(["2024-01-01", "2024-01-02"], type=f.type)
        elif f.name == "qpx_version":
            arrays[f.name] = pa.array(["1.0", "1.0"], type=f.type)
        else:
            arrays[f.name] = pa.nulls(2, type=f.type)
    table = pa.table(arrays, schema=ds_schema)
    result = DatasetSchema.validate_full(table)
    pk_issues = [i for i in result.issues if i.check == "duplicate_pk"]
    assert len(pk_issues) == 1
    assert pk_issues[0].severity == "warning"


def test_validate_backward_compat():
    """validate() returns list[str] for backward compatibility."""
    schema = FeatureSchema.get_arrow_schema()

    # With missing column
    fields_to_keep = [f for f in schema if f.name != "sequence"]
    arrays = {f.name: pa.nulls(1, type=f.type) for f in fields_to_keep}
    table = pa.table(arrays, schema=pa.schema(fields_to_keep))
    errors = FeatureSchema.validate(table)
    assert isinstance(errors, list)
    assert all(isinstance(e, str) for e in errors)
    assert any("sequence" in e for e in errors)

    # Valid table
    arrays = {f.name: pa.nulls(1, type=f.type) for f in schema}
    table = pa.table(arrays, schema=schema)
    assert FeatureSchema.validate(table) == []


def test_dataset_validate(dataset_dir):
    """Dataset.validate: all structures, specific structure, missing structure."""
    import qpx

    with qpx.open_dataset(dataset_dir) as ds:
        # All structures
        results = ds.validate()
        assert len(results) > 0
        for name, result in results.items():
            assert isinstance(result, ValidationResult)
            assert result.is_valid, f"{name}: {result.summary}"

        # Specific structure
        results = ds.validate(structures=["feature"])
        assert "feature" in results
        assert len(results) == 1
        assert results["feature"].is_valid

        # Missing structure
        results = ds.validate(structures=["mz"])
        assert "mz" in results
        assert not results["mz"].is_valid
        assert any(i.check == "missing_structure" for i in results["mz"].issues)


def test_cli_validate(dataset_dir, feature_parquet):
    """CLI validate: dataset, specific structure, single file, no args error."""
    from qpx.cli.main import qpx_main

    runner = CliRunner()

    result = runner.invoke(qpx_main, ["validate", "--dataset-path", str(dataset_dir)])
    assert result.exit_code == 0
    assert "VALID" in result.output

    result = runner.invoke(
        qpx_main,
        ["validate", "--dataset-path", str(dataset_dir), "--structure", "feature"],
    )
    assert result.exit_code == 0
    assert "feature" in result.output

    result = runner.invoke(qpx_main, ["validate", "--file", str(feature_parquet)])
    assert result.exit_code == 0
    assert "VALID" in result.output

    result = runner.invoke(qpx_main, ["validate"])
    assert result.exit_code != 0
