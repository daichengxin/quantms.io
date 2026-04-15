# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2025-03-22

### Added

- **CLI command groups**: `convert`, `transform`, `query`, `info`, `validate`, `ontology`
- **Converters**: DIA-NN, MaxQuant, QuantMS (mzTab LFQ & TMT), FragPipe, mzIdentML, SDRF
- **Transform commands**:
  - `gene-map` — map gene names from FASTA to QPX parquet data
  - `quantify` — protein-level quantification via mokume (DirectLFQ, MaxLFQ, iBAQ, TopN, Sum)
  - `normalize-accessions` — normalize protein accession formats (full ↔ bare UniProt)
  - `update-metadata` — update sample/run metadata from a revised SDRF
- **Query commands**: `sql`, `filter`, `head` for interactive dataset exploration
- **Info commands**: dataset summary, Arrow schema display, Parquet metadata inspection
- **Validate command**: schema validation with column presence, type matching, null checks, and PK uniqueness
- **Ontology management**: `info`, `update`, `build`, `search` for PSI-MS and PRIDE CV terms
- **Python API**: `qpx.open_dataset()`, `qpx.read_feature()`, `qpx.read_psm()`, `qpx.read_pg()`, etc.
- **Dataset class**: unified access to all QPX structures with DuckDB-backed SQL queries
- **DatasetCollection**: work with multiple datasets simultaneously
- **Writers**: Parquet writers for all QPX structures (Feature, PSM, PG, Sample, Run, MzSpectra, PepMap, Ontology, Provenance)
- **Views**: analytical views for protein, peptide, QC, and sample summaries
- **Schema validation**: YAML-defined canonical schemas for all data structures
- **Score normalization**: automatic PSI-MS ontology mapping for search engine scores
- **PRIDE API integration**: `--enrich-pride` flag for automatic project metadata enrichment
- **S3 support**: read QPX datasets from S3-compatible storage
- **MkDocs documentation**: full CLI reference with auto-generated parameter tables and usage examples

### Changed

- Replaced per-converter `_safe_float_val` with shared `safe_float` utility
- Optimized `_table_exists` to use parameterized `information_schema` query
- Increased DIA-NN `file_num` default from 50 to 100 for larger cohort support
- Defensive handling of None/NaN in gene mapping (`_resolve_gene_names`, `annotate_dataframe`)

### Fixed

- `.gitignore` duplicate entry removal
- SQL injection safety in `_table_exists` (parameterized query)
- Gene mapping crash on None/NaN protein accessions
- Gene mapping `list(dict)` producing keys instead of preserving struct

## [0.0.4] - 2025-03-01

### Added

- Initial PyPI release
- Basic converter framework for QuantMS, DIA-NN, MaxQuant
- QPX Parquet schema definition
- CLI entry point (`qpxc`)
