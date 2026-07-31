# Sample Metadata

`sample.parquet` stores one row per biological sample, capturing the biological properties of each specimen in the study. It is derived from the community [SDRF](https://github.com/bigbio/proteomics-sample-metadata) (Sample and Data Relationship Format) standard, with bracket-style column names translated to clean `snake_case` field names during ingestion.

See the full YAML schema in [`sample.yaml`](schemas/sample.yaml).

!!! info "SDRF as the ingestion format"
    QPX does **not** invent its own sample metadata model. It adopts the community SDRF standard and converts it to query-friendly Parquet files. At ingestion time, bracket-style column names are translated to clean `snake_case` field names and the data is split by entity type.

---

## PyArrow Schema

```python
import pyarrow as pa

sample_schema = pa.schema(
    [
        # Identity -- maps to SDRF "source name"
        pa.field("sample_accession", pa.string()),  # PK, required
        # Mandatory biological properties
        pa.field("organism", pa.string()),  # required
        pa.field("organism_part", pa.string()),  # required
        # Optional standard properties
        pa.field("disease", pa.string(), nullable=True),
        pa.field("cell_line", pa.string(), nullable=True),
        pa.field("cell_type", pa.string(), nullable=True),
        pa.field("sex", pa.string(), nullable=True),
        pa.field("age", pa.string(), nullable=True),
        pa.field("developmental_stage", pa.string(), nullable=True),
        pa.field("ancestry", pa.string(), nullable=True),
        pa.field("individual", pa.string(), nullable=True),
        pa.field("sample_description", pa.string(), nullable=True),
        # Extra columns from SDRF characteristics are added dynamically
        # e.g. pa.field("biological_replicate", pa.string(), nullable=True),
        # e.g. pa.field("bmi", pa.string(), nullable=True),
    ]
)
```

## Field Reference

### Core fields (always present)

| Field | Description | Type | Required |
|---|---|---|---|
| `sample_accession` | Unique sample identifier (= SDRF `source name`) | string | yes |
| `organism` | Species (e.g. `"Homo sapiens"`) | string | yes |
| `organism_part` | Tissue or organ (e.g. `"Brain; Frontal cortex"`) | string | yes |
| `disease` | Disease state (e.g. `"breast cancer"`) | string | no |
| `cell_line` | Cell line name (e.g. `"HeLa"`) | string | no |
| `cell_type` | Cell type (e.g. `"T cell"`) | string | no |
| `sex` | Biological sex (e.g. `"female"`) | string | no |
| `age` | Age of the specimen (e.g. `"45"`) | string | no |
| `developmental_stage` | Developmental stage | string | no |
| `ancestry` | Ancestry category | string | no |
| `individual` | Individual/patient identifier | string | no |
| `sample_description` | Free-text sample description | string | no |

### Extra columns (dataset-specific)

Any non-standard SDRF `characteristics[X]` column becomes its own string column in the Parquet file. For example:

| Extra column | SDRF source |
|---|---|
| `biological_replicate` | `characteristics[biological replicate]` |
| `bmi` | `characteristics[bmi]` |
| `treatment` | `characteristics[treatment]` |

!!! note "Multi-valued fields"
    When a sample has multiple values for a property (e.g., a spike-in with two organisms), the values are joined with `"; "` into a single string: `"Homo sapiens; Saccharomyces cerevisiae"`.

## Example Rows

| sample_accession | organism | organism_part | disease | sex | biological_replicate |
|---|---|---|---|---|---|
| Sample_01 | Homo sapiens | breast | breast cancer | female | 1 |
| Sample_02 | Homo sapiens | breast | normal | female | 2 |
| Sample_03 | Homo sapiens; Saccharomyces cerevisiae | breast | breast cancer | male | 3 |

!!! tip "Reading sample.parquet"
    ```python
    import pyarrow.parquet as pq

    samples = pq.read_table("PXD014414.sample.parquet")
    print(samples.column_names)
    # ['sample_accession', 'organism', 'organism_part', 'disease', ..., 'biological_replicate']

    # Get unique organisms
    import pyarrow.compute as pc

    organisms = pc.unique(samples.column("organism")).to_pylist()
    print(organisms)  # ['Homo sapiens', 'Homo sapiens; Saccharomyces cerevisiae']
    ```

!!! tip "Querying with DuckDB"
    ```sql
    -- All unique organisms in the project
    SELECT DISTINCT organism
    FROM 'PXD014414.sample.parquet';

    -- Samples with a specific disease
    SELECT sample_accession, organism_part, disease
    FROM 'PXD014414.sample.parquet'
    WHERE disease LIKE '%breast cancer%';

    -- Group by extra column
    SELECT biological_replicate, COUNT(*) AS n_samples
    FROM 'PXD014414.sample.parquet'
    GROUP BY biological_replicate;
    ```

---

## Related Pages

- [Run Metadata](run.md) -- `run.parquet`
- [Dataset Metadata](dataset.md) -- `dataset.parquet`
- [SDRF Conversion](sdrf-mapping.md) -- Column mapping and ingestion process
- [QPX Format Overview](index.md)
