# Versioning

QPX uses a simple versioning scheme to track specification changes and ensure that tools can correctly interpret files produced by different versions of the format.

## Version Format

QPX versions follow a **`{major}.{minor}`** format -- for example, **1.0**.

This is a two-component scheme without patch numbers. The version applies to the entire QPX specification, not to individual views.

## Rules

- **Major releases** (`X.0`): Introduce backward-incompatible changes. These may include removing fields, changing field types, renaming views, or restructuring the data model. Tools built for a previous major version may not be able to read files produced by a new major version without updates.
- **Minor updates** (`1.X`): Introduce backward-compatible additions only. These may include adding new optional fields, new views, or new metadata keys. Tools built for an earlier minor version within the same major version can safely read files produced by a newer minor version (unknown fields are ignored).

!!! note
    QPX does not use patch versions. Bug fixes to the specification text that do not change the data model are tracked in the documentation changelog but do not increment the version number.

## Version in Files

All QPX views (PSM, Feature, PG, Peptide, Protein, MZ, and others serialized as Parquet) include a `qpx_version` field in their file-level metadata. This metadata is stored as a key-value pair in the Parquet file footer and identifies which version of the QPX specification was used to generate the file.

```python
import pyarrow.parquet as pq

# Writing a file with version metadata
metadata = {'qpx_version': '1.0'}
```

!!! tip "Reading the version from a file"
    ```python
    import pyarrow.parquet as pq

    parquet_file = pq.ParquetFile("experiment.psm.parquet")
    schema_metadata = parquet_file.schema_arrow.metadata
    qpx_version = schema_metadata.get(b'qpx_version', b'unknown').decode()
    print(f"QPX version: {qpx_version}")
    ```

## Current Version

The current QPX specification version is **1.0**.

This version defines all core views (PSM, Feature, PG, Peptide, Protein, MZ), expression views (Absolute, Differential), and metadata views (SDRF, Project).

## Software Provider

In addition to the specification version, every QPX file records the software that generated the data. This is stored in the `software_provider` metadata field in the Parquet file footer. The value identifies the tool name and version that produced the file:

```json
{
  "software_provider": {
    "name": "quantms",
    "version": "1.3.0"
  }
}
```

!!! info "Why track the software provider?"
    Tracking the generating software enables reproducibility and debugging. If a downstream analysis produces unexpected results, the software provider metadata allows users to determine which tool version created the input data and whether known issues apply.

The `software_provider` field is a free-text string in the Parquet metadata. The JSON structure shown above is the recommended format, but tools may also use a flat string such as `"QuantMS 1.3.0"`.

## Version Compatibility Matrix

| Reader Version | File Version | Compatible? | Notes |
|---------------|-------------|-------------|-------|
| 1.0 | 1.0 | Yes | Exact match |
| 1.1 | 1.0 | Yes | Reader is newer, can handle older files |
| 1.0 | 1.1 | Yes | Reader ignores unknown optional fields |
| 2.0 | 1.0 | Maybe | Major version change; check migration guide |
| 1.0 | 2.0 | No | Reader cannot handle breaking changes |

## See Also

- [Serialization](serialization.md) -- serialization formats and file-level metadata details
- [QPX Format Overview](index.md) -- overview of all QPX views and the current specification status
