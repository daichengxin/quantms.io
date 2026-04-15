"""Tests for QPX schema layer: FieldDef, ViewSchema, and all schema classes."""

import pyarrow as pa

from qpx.core.data import (
    DatasetSchema,
    FeatureSchema,
    FieldDef,
    MzSchema,
    OntologySchema,
    PgSchema,
    ProvenanceSchema,
    PsmSchema,
    RunSchema,
    SampleSchema,
    ViewSchema,
)

ALL_SCHEMAS = [
    FeatureSchema,
    PsmSchema,
    PgSchema,
    MzSchema,
    SampleSchema,
    RunSchema,
    DatasetSchema,
    OntologySchema,
    ProvenanceSchema,
]


def test_field_def():
    """FieldDef: nullable, optional, types, doc metadata, complex nested types."""
    # Required field not nullable
    f = FieldDef("my_col", pa.string(), nullable=False, doc="test")
    arrow = f.to_arrow_field()
    assert arrow.name == "my_col"
    assert arrow.type == pa.string()
    assert arrow.nullable is False

    # Nullable field
    f = FieldDef("score", pa.float64(), nullable=True, doc="test nullable")
    assert f.to_arrow_field().nullable is True

    # Optional field
    f = FieldDef("opt_col", pa.int32(), nullable=True, optional=True, doc="optional")
    assert f.to_arrow_field().nullable is True
    assert f.optional is True

    # Default field is nullable
    f = FieldDef("default_col", pa.string())
    assert f.to_arrow_field().nullable is True

    # Preserves arrow type
    f = FieldDef("scan", pa.list_(pa.int32()), nullable=False)
    assert f.to_arrow_field().type == pa.list_(pa.int32())

    # Doc in metadata
    f = FieldDef("documented_col", pa.string(), nullable=False, doc="My documentation")
    arrow = f.to_arrow_field()
    assert arrow.metadata is not None
    assert arrow.metadata[b"doc"] == b"My documentation"

    # No doc, no metadata
    f = FieldDef("no_doc_col", pa.string(), nullable=False)
    assert f.to_arrow_field().metadata is None

    # Complex nested type
    score_type = pa.list_(
        pa.struct(
            [
                pa.field("score_name", pa.string()),
                pa.field("score_value", pa.float64()),
                pa.field("higher_better", pa.bool_(), nullable=True),
            ]
        )
    )
    f = FieldDef("scores", score_type, nullable=True)
    assert f.to_arrow_field().type == score_type


def test_view_schema_basic():
    """ViewSchema: arrow schema, iter_fields, field names match, is_optional."""
    schema = FeatureSchema.get_arrow_schema()
    assert isinstance(schema, pa.Schema)

    fields = list(FeatureSchema._iter_fields())
    assert len(fields) > 0
    for name, field in fields:
        assert isinstance(name, str)
        assert isinstance(field, FieldDef)

    # Field names match
    field_names_from_schema = set(schema.names)
    field_names_from_iter = {name for name, _ in FeatureSchema._iter_fields()}
    assert field_names_from_schema == field_names_from_iter

    # is_optional
    assert FeatureSchema._is_optional("pg_global_qvalue") is True
    assert FeatureSchema._is_optional("gg_accessions") is True
    assert FeatureSchema._is_optional("sequence") is False
    assert FeatureSchema._is_optional("nonexistent_field") is False


def test_schema_validation():
    """Validation: valid table, missing required, type mismatch, optional absent."""
    schema = FeatureSchema.get_arrow_schema()

    # Valid table returns no errors
    arrays = {f.name: pa.nulls(1, type=f.type) for f in schema}
    table = pa.table(arrays, schema=schema)
    assert FeatureSchema.validate(table) == []

    # Missing required column
    fields_to_keep = [f for f in schema if f.name != "sequence"]
    arrays = {f.name: pa.nulls(1, type=f.type) for f in fields_to_keep}
    table = pa.table(arrays, schema=pa.schema(fields_to_keep))
    errors = FeatureSchema.validate(table)
    assert len(errors) >= 1
    assert any("sequence" in e for e in errors)

    # Type mismatch
    wrong_fields = [pa.field(f.name, pa.string()) if f.name == "charge" else f for f in schema]
    arrays = {f.name: pa.nulls(1, type=f.type) for f in wrong_fields}
    arrays["charge"] = pa.array(["wrong_type"], type=pa.string())
    table = pa.table(arrays, schema=pa.schema(wrong_fields))
    errors = FeatureSchema.validate(table)
    assert any("charge" in e and "int16" in e for e in errors)

    # Optional column absent is fine
    fields_to_keep = [f for f in schema if f.name != "pg_global_qvalue"]
    arrays = {f.name: pa.nulls(1, type=f.type) for f in fields_to_keep}
    table = pa.table(arrays, schema=pa.schema(fields_to_keep))
    errors = FeatureSchema.validate(table)
    assert not any("pg_global_qvalue" in e for e in errors)


def test_all_schemas_have_valid_primary_keys():
    """Every schema has a _primary_key tuple whose fields exist in the schema."""
    for schema_inst in ALL_SCHEMAS:
        assert isinstance(schema_inst._primary_key, tuple)
        assert len(schema_inst._primary_key) > 0
        schema = schema_inst.get_arrow_schema()
        for pk_field in schema_inst._primary_key:
            assert pk_field in schema.names, f"PK field '{pk_field}' not in {schema_inst.__name__}"


def test_yaml_loader_and_nested_types():
    """YAML loader returns ViewSchema; nested structs resolve; sample extension works."""
    from qpx.core.data.loader import load_schema

    schema = load_schema("feature")
    assert isinstance(schema, ViewSchema)

    # Nested struct: modifications = list<struct<name, accession, positions>>
    arrow_schema = FeatureSchema.get_arrow_schema()
    mod_type = arrow_schema.field("modifications").type
    assert isinstance(mod_type, pa.ListType)
    inner = mod_type.value_type
    assert isinstance(inner, pa.StructType)
    assert inner.get_field_index("name") >= 0
    assert inner.get_field_index("positions") >= 0

    # Nested struct: intensities = list<struct<label, intensity>>
    int_type = arrow_schema.field("intensities").type
    assert isinstance(int_type, pa.ListType)
    inner = int_type.value_type
    assert isinstance(inner, pa.StructType)
    assert inner.get_field_index("label") >= 0
    assert inner.get_field_index("intensity") >= 0

    # Sample extra columns
    extended = SampleSchema.get_extended_arrow_schema(["bmi", "biological_replicate"])
    assert "bmi" in extended.names
    assert extended.field("bmi").type == pa.string()
    assert extended.field("organism").type == pa.string()
