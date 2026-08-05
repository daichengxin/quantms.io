"""Regression tests for identity-preserving gene annotation writes."""

import pyarrow.parquet as pq

from qpx import Dataset
from qpx.transforms.gene_mapping import GeneMappingTransform
from qpx.writers.feature import FeatureWriter
from tests.conftest import make_feature_record


def test_write_annotated_features_preserves_source_identity(tmp_path, monkeypatch):
    """Annotation must not replace a producer-specific Feature identity recipe."""
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    source_path = source_dir / "fragpipe.feature.parquet"
    records = []
    for voltage in (-65.0, -45.0):
        record = make_feature_record(run_file_name="experiment", scan=[])
        record.update(
            {
                "rt": None,
                "quantification_unit_id": "experiment",
                "compensation_voltage": voltage,
            }
        )
        records.append(record)

    composite = ("quantification_unit_id", "peptidoform", "charge", "compensation_voltage")
    with FeatureWriter(source_path, identity_composite=composite) as writer:
        writer.write_batch(records)

    dataset = Dataset(source_dir, structures=["feature"])
    fasta = tmp_path / "empty.fasta"
    fasta.touch()
    transform = GeneMappingTransform(fasta, fetch_accessions=False)
    source_frame = dataset.feature.to_df()
    monkeypatch.setattr(transform, "annotate_dataset_features", lambda *_args, **_kwargs: source_frame)

    output_path = tmp_path / "annotated.feature.parquet"
    transform.write_annotated_features(dataset, output_path)

    source = pq.read_table(source_path)
    output = pq.read_table(output_path)
    assert output.column("feature_id").to_pylist() == source.column("feature_id").to_pylist()
    assert output.schema.metadata[b"identity_composite"] == b",".join(field.encode() for field in composite)
