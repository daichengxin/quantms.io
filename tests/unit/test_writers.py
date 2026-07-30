"""Tests for QPX writer layer: BaseWriter subclasses and round-trip validation."""

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from qpx.core.data import FeatureSchema, PsmSchema
from qpx.core.parquet_io import parquet_row_count, read_parquet_metadata
from qpx.writers import (
    DatasetWriter,
    FeatureWriter,
    OntologyWriter,
    PgWriter,
    ProvenanceWriter,
    PsmWriter,
    RunWriter,
    SampleWriter,
)
from tests.conftest import (
    make_dataset_record,
    make_feature_record,
    make_ontology_record,
    make_pg_record,
    make_provenance_record,
    make_psm_record,
    make_run_record,
    make_sample_record,
)


def test_writers_create_valid_parquet(tmp_path):
    """All writers produce readable Parquet files with correct row count."""
    configs = [
        (FeatureWriter, make_feature_record, "feature"),
        (PsmWriter, make_psm_record, "psm"),
        (PgWriter, make_pg_record, "pg"),
        (SampleWriter, make_sample_record, "sample"),
        (RunWriter, make_run_record, "run"),
        (DatasetWriter, make_dataset_record, "dataset"),
        (OntologyWriter, make_ontology_record, "ontology"),
        (ProvenanceWriter, make_provenance_record, "provenance"),
    ]
    for writer_cls, record_fn, name in configs:
        path = tmp_path / f"test.{name}.parquet"
        with writer_cls(path) as w:
            w.write_batch([record_fn()])
        assert path.exists()
        table = pq.read_table(path)
        assert table.num_rows == 1


def test_writer_footer_metadata(tmp_path):
    """Footer metadata contains file_type, version, creator, software_provider, scan_format."""
    path = tmp_path / "test.feature.parquet"
    with FeatureWriter(path, creator="test_suite", software_provider="my_tool", scan_format="native") as w:
        w.write_batch([make_feature_record()])

    meta = read_parquet_metadata(path)
    assert meta["file_type"] == "feature_file"
    assert meta["creator"] == "test_suite"
    assert "qpx_version" in meta
    assert "uuid" in meta
    assert "creation_date" in meta
    assert meta["software_provider"] == "my_tool"
    assert meta["scan_format"] == "native"

    # Metadata should be in footer, not in data columns
    table = pq.read_table(path)
    data_columns = table.schema.names
    assert "qpx_version" not in data_columns
    assert "file_type" not in data_columns


def test_writer_schema_validation(tmp_path):
    """Writer rejects bad types and accepts valid tables."""
    # Reject bad types
    path = tmp_path / "bad.feature.parquet"
    writer = FeatureWriter(path)
    try:
        schema = FeatureSchema.get_arrow_schema()
        wrong_fields = [pa.field("charge", pa.string()) if f.name == "charge" else f for f in schema]
        wrong_schema = pa.schema(wrong_fields)
        arrays = {f.name: pa.nulls(1, type=f.type) for f in wrong_schema}
        arrays["charge"] = pa.array(["wrong"], type=pa.string())
        table = pa.table(arrays, schema=wrong_schema)
        with pytest.raises(ValueError, match="Schema validation failed"):
            writer.write_table(table)
    finally:
        writer.close()

    # Accept valid table
    path2 = tmp_path / "good.feature.parquet"
    writer = FeatureWriter(path2)
    schema = writer.arrow_schema
    batch = pa.RecordBatch.from_pylist([make_feature_record()], schema=schema)
    table = pa.Table.from_batches([batch])
    with writer:
        writer.write_table(table)
    assert parquet_row_count(path2) == 1


def test_writers_support_de_novo_records_without_database_fields(tmp_path):
    """PSM and feature writers accept records from a database-free workflow."""
    psm_record = make_psm_record(sequence="DENOVO", peptidoform="DENOVO")
    psm_record["protein_accessions"] = None

    psm_path = tmp_path / "denovo.psm.parquet"
    with PsmWriter(psm_path) as writer:
        writer.write_batch([psm_record])
    psm_table = pq.read_table(psm_path)
    assert psm_table.column("is_decoy").to_pylist() == [False]
    assert psm_table.column("protein_accessions").to_pylist() == [None]
    assert PsmSchema.validate(psm_table) == []

    feature_record = make_feature_record(sequence="DENOVO", peptidoform="DENOVO")
    feature_record["anchor_protein"] = None
    feature_record["pg_accessions"] = None
    feature_record["unique"] = None

    feature_path = tmp_path / "denovo.feature.parquet"
    with FeatureWriter(feature_path) as writer:
        writer.write_batch([feature_record])
    feature_table = pq.read_table(feature_path)
    assert feature_table.column("is_decoy").to_pylist() == [False]
    assert feature_table.column("anchor_protein").to_pylist() == [None]
    assert FeatureSchema.validate(feature_table) == []


def test_writer_batching(tmp_path):
    """Batch flush, remaining buffer flush on close, and multiple write_batch calls."""
    # Batch size flush
    path = tmp_path / "batch.feature.parquet"
    with FeatureWriter(path, batch_size=2) as w:
        records = [make_feature_record(sequence=f"SEQ{i}") for i in range(5)]
        w.write_batch(records)
    assert parquet_row_count(path) == 5

    # Remaining buffer flushed on close
    path2 = tmp_path / "buffer.feature.parquet"
    with FeatureWriter(path2, batch_size=100) as w:
        records = [make_feature_record(sequence=f"SEQ{i}") for i in range(3)]
        w.write_batch(records)
    assert parquet_row_count(path2) == 3

    # Multiple write_batch calls
    path3 = tmp_path / "multi.feature.parquet"
    with FeatureWriter(path3) as w:
        w.write_batch([make_feature_record(sequence="SEQ1")])
        w.write_batch([make_feature_record(sequence="SEQ2")])
        w.write_batch([make_feature_record(sequence="SEQ3")])
    assert parquet_row_count(path3) == 3


def test_writer_compression(tmp_path):
    """Default compression is zstd; custom compression is recorded."""
    path = tmp_path / "default.feature.parquet"
    with FeatureWriter(path) as w:
        w.write_batch([make_feature_record()])
    assert read_parquet_metadata(path)["compression_format"] == "zstd"

    path2 = tmp_path / "snappy.feature.parquet"
    with FeatureWriter(path2, compression="snappy") as w:
        w.write_batch([make_feature_record()])
    assert read_parquet_metadata(path2)["compression_format"] == "snappy"
