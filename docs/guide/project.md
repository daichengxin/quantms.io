# Project Metadata Management

!!! warning "v2.0 Migration"
    In QPX v2.0, project metadata is managed through the Dataset API. The previous command-line tools for project management have been replaced with a simplified API-based approach. Project-level metadata is now stored in [`dataset.parquet`](../spec/dataset.md), [`sample.parquet`](../spec/sample.md), and [`run.parquet`](../spec/run.md).

## Overview

QPX v2.0 manages project metadata through the Dataset API, providing a unified interface for accessing and managing dataset-level information, sample metadata, and run-level details.

## Accessing Project Metadata

Load a dataset and access project-level metadata:

```python
import qpx

# Load dataset
ds = qpx.Dataset("path/to/dataset/")

# Access project-level metadata
print(ds.dataset_meta)

# Access specific metadata fields
if hasattr(ds.dataset_meta, "name"):
    print(f"Dataset name: {ds.dataset_meta.name}")
if hasattr(ds.dataset_meta, "description"):
    print(f"Description: {ds.dataset_meta.description}")
if hasattr(ds.dataset_meta, "version"):
    print(f"Version: {ds.dataset_meta.version}")
```

## Dataset Metadata Structure

The `dataset_meta` attribute provides access to dataset-level information:

```python
import qpx

ds = qpx.Dataset("./output/")

# Dataset metadata includes:
# - name: Dataset name or identifier
# - description: Dataset description
# - version: QPX format version
# - created_at: Creation timestamp
# - modified_at: Last modification timestamp
# - software: Software used to generate the data
# - parameters: Processing parameters

# Access as attributes
metadata = ds.dataset_meta
print(f"Name: {metadata.name}")
print(f"Version: {metadata.version}")
```

## Sample Metadata

Access sample-level metadata through the Dataset API:

```python
import qpx

ds = qpx.Dataset("./output/")

# Sample metadata is stored in sample.parquet
if hasattr(ds, "sample"):
    sample_data = ds.sample.data
    print(sample_data.head())

    # Access specific sample information
    sample_ids = sample_data["sample_accession"].unique()
    print(f"Number of samples: {len(sample_ids)}")

    # Sample metadata includes:
    # - sample_accession: Unique sample identifier
    # - organism: Organism(s) studied
    # - organism_part: Tissue or cell type
    # - disease: Disease state (if applicable)
    # - biological_replicate: Replicate number
    # - condition: Experimental condition
```

## Run Metadata

Access MS run-level metadata:

```python
import qpx

ds = qpx.Dataset("./output/")

# Run metadata is stored in run.parquet
if hasattr(ds, "run"):
    run_data = ds.run.data
    print(run_data.head())

    # Access specific run information
    run_ids = run_data["run_accession"].unique()
    print(f"Number of runs: {len(run_ids)}")

    # Run metadata includes:
    # - run_accession: Unique run identifier
    # - sample_accession: Associated sample
    # - instrument: Mass spectrometer used
    # - acquisition_date: When the run was acquired
    # - raw_file_name: Original raw file name
```

## Viewing Complete Metadata

Generate a comprehensive metadata summary:

```python
import qpx


def print_dataset_summary(dataset_path):
    """Print comprehensive dataset metadata summary."""
    ds = qpx.Dataset(dataset_path)

    print("=" * 60)
    print("QPX Dataset Metadata Summary")
    print("=" * 60)
    print()

    # Dataset-level metadata
    print("Dataset Information:")
    if hasattr(ds, "dataset_meta"):
        if hasattr(ds.dataset_meta, "name"):
            print(f"  Name: {ds.dataset_meta.name}")
        if hasattr(ds.dataset_meta, "version"):
            print(f"  Version: {ds.dataset_meta.version}")
        if hasattr(ds.dataset_meta, "description"):
            print(f"  Description: {ds.dataset_meta.description}")
    print()

    # Sample information
    if hasattr(ds, "sample") and ds.sample.count() > 0:
        print("Sample Information:")
        print(f"  Total samples: {ds.sample.count()}")
        if "organism" in ds.sample.data.columns:
            organisms = ds.sample.data["organism"].unique()
            print(f"  Organisms: {', '.join(organisms)}")
        if "condition" in ds.sample.data.columns:
            conditions = ds.sample.data["condition"].unique()
            print(f"  Conditions: {', '.join(map(str, conditions))}")
    print()

    # Run information
    if hasattr(ds, "run") and ds.run.count() > 0:
        print("Run Information:")
        print(f"  Total runs: {ds.run.count()}")
        if "instrument" in ds.run.data.columns:
            instruments = ds.run.data["instrument"].unique()
            print(f"  Instruments: {', '.join(instruments)}")
    print()

    # Data availability
    print("Data Availability:")
    print(f"  PSM data: {'Yes' if hasattr(ds, 'psm') and ds.psm.count() > 0 else 'No'}")
    print(f"  Feature data: {'Yes' if hasattr(ds, 'feature') and ds.feature.count() > 0 else 'No'}")
    print(f"  Protein group data: {'Yes' if hasattr(ds, 'pg') and ds.pg.count() > 0 else 'No'}")
    print()

    print("=" * 60)


# Usage
print_dataset_summary("./output/")
```

## Metadata from SDRF Files

QPX can parse sample metadata from SDRF (Sample and Data Relationship Format) files:

```python
import qpx
from qpx.core.sdrf import SDRFHandler

# Parse SDRF file
sdrf = SDRFHandler("./metadata.sdrf.tsv")

# Access SDRF metadata
print("SDRF Metadata:")
print(f"  Samples: {len(sdrf.samples)}")
print(f"  Runs: {len(sdrf.runs)}")

# Get sample information
for sample in sdrf.samples:
    print(f"  Sample: {sample.sample_accession}")
    print(f"    Organism: {sample.organism}")
    print(f"    Condition: {sample.condition}")
```

## Provenance Information

Track data processing provenance:

```python
import qpx

ds = qpx.Dataset("./output/")

# Provenance information is stored in provenance.parquet
if hasattr(ds, "provenance"):
    prov_data = ds.provenance.data
    print("Processing Provenance:")
    print(prov_data)

    # Provenance includes:
    # - software_name: Software used
    # - software_version: Version number
    # - parameters: Processing parameters
    # - timestamp: When processing occurred
    # - input_files: Source data files
```

## Integration with PRIDE Archive

For datasets from PRIDE Archive, metadata can be enriched with public repository information. While the command-line tools have been removed, you can use the PRIDE MCP (Model Context Protocol) tools to fetch metadata:

```python
# Example of accessing PRIDE metadata
# Note: This requires the PRIDE MCP server to be configured

# Fetch project details
# project_details = get_pride_project_details("PXD001234")

# The returned metadata includes:
# - title: Project title
# - description: Project description
# - organism: Studied organisms
# - instrument: Mass spectrometers used
# - publication: Associated publication
# - submission_date: When submitted to PRIDE

# This metadata can be used to enrich your QPX dataset
```

## Best Practices

### Metadata Documentation

Keep metadata well-documented:

```python
import qpx

# When creating datasets, include comprehensive metadata
# This helps with reproducibility and data sharing

ds = qpx.Dataset("./output/")

# Document:
# 1. Data source (instrument, date, operator)
# 2. Processing software and versions
# 3. Search parameters
# 4. Sample preparation details
# 5. Experimental design
```

### Metadata Validation

Validate metadata completeness:

```python
import qpx


def validate_metadata(dataset_path):
    """Check for required metadata fields."""
    ds = qpx.Dataset(dataset_path)

    issues = []

    # Check dataset-level metadata
    if not hasattr(ds, "dataset_meta"):
        issues.append("Missing dataset metadata")
    elif not hasattr(ds.dataset_meta, "name"):
        issues.append("Missing dataset name")

    # Check sample metadata
    if not hasattr(ds, "sample") or ds.sample.count() == 0:
        issues.append("Missing sample metadata")

    # Check run metadata
    if not hasattr(ds, "run") or ds.run.count() == 0:
        issues.append("Missing run metadata")

    # Report
    if issues:
        print("Metadata Validation Issues:")
        for issue in issues:
            print(f"  - {issue}")
        return False
    else:
        print("Metadata validation passed!")
        return True


# Usage
validate_metadata("./output/")
```

### Metadata Export

Export metadata for sharing:

```python
import qpx
import json


def export_metadata(dataset_path, output_file):
    """Export dataset metadata to JSON."""
    ds = qpx.Dataset(dataset_path)

    metadata = {}

    # Dataset-level metadata
    if hasattr(ds, "dataset_meta"):
        metadata["dataset"] = {
            "name": getattr(ds.dataset_meta, "name", None),
            "version": getattr(ds.dataset_meta, "version", None),
            "description": getattr(ds.dataset_meta, "description", None),
        }

    # Sample summary
    if hasattr(ds, "sample") and ds.sample.count() > 0:
        metadata["samples"] = {"count": ds.sample.count(), "sample_ids": ds.sample.data["sample_accession"].tolist()}

    # Run summary
    if hasattr(ds, "run") and ds.run.count() > 0:
        metadata["runs"] = {"count": ds.run.count(), "run_ids": ds.run.data["run_accession"].tolist()}

    # Save to file
    with open(output_file, "w") as f:
        json.dump(metadata, f, indent=2)

    print(f"Metadata exported to: {output_file}")


# Usage
export_metadata("./output/", "./metadata_export.json")
```

## Migration from v1.0

If you have project.json files from QPX v1.0, the metadata is now distributed across the new format:

| v1.0 (project.json)      | v2.0 Location          |
| ------------------------ | ---------------------- |
| `accession`              | `dataset_meta.name`    |
| `title`                  | `dataset_meta.name`    |
| `description`            | `dataset_meta.description` |
| `samples[]`              | `sample.parquet`       |
| `organism`               | `sample.parquet`       |
| `instrument`             | `run.parquet`          |
| `software`               | `provenance.parquet`   |
| `quantification_method`  | `dataset_meta`         |
| `publication`            | `dataset_meta`         |

## Related Documentation

- [Dataset Specification](../spec/dataset.md) - Dataset metadata format
- [Sample Specification](../spec/sample.md) - Sample metadata format
- [Run Specification](../spec/run.md) - Run metadata format
- [Provenance Specification](../spec/provenance.md) - Processing provenance
- [SDRF Mapping](../spec/sdrf-mapping.md) - SDRF to QPX mapping
