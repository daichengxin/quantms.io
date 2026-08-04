# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Spectronaut converter**: `qpxc convert spectronaut` — full support for Spectronaut report TSV files, producing feature.parquet and pg.parquet with DuckDB-accelerated batch processing
- **CPTAC CDAP converter**: `qpxc convert cdap` — convert CPTAC CDAP `.psm` study directories to QPX psm/feature/pg/dataset/ontology/provenance views
- **Full-spectra mz converter**: `qpxc convert mz` — convert a directory of mzML / `.mzML.gz` files to a single `mz.parquet`; each spectrum carries `run_file_name` + `scan` for linking back to PSM/feature
- **pdc2qpx pipeline**: `qpxc pdc2qpx` — one-shot PDC/CPTAC download (via pridepy, `qpx[pdc]` extra) + CDAP + full-spectra conversion into an entire QPX dataset
- **Shared channel-label resolution**: `qpx/converters/channel_labels.py` — single source for canonical TMT/iTRAQ/LFQ labels via sdrf-pipelines `channel_map`, used by both QuantMS mzTab and OpenMS `-out_qpx` paths
- **openms-consensus interim protein intensity**: the `openms-consensus` converter now fills `pg.intensity` with an interim, **unnormalized sum of each group's unique peptides** per `(protein group, grouped_runs, label)` (the quantms `unique_peptides` policy) instead of leaving it null, until OpenMS `-out_qpx` ships the authoritative quant. Every quantified row is stamped with a `quantification_method` cv_param; `--pg-top N` bounds the peptides used (`0` = all; `3` mirrors the ProteomicsLFQ/IsobaricWorkflow default)
- **openms-consensus channel/SDRF consistency check**: when an SDRF is provided, the converter compares the isobaric channels read from the consensusXML maps against the SDRF `comment[label]` set and logs a warning for any channel present in one but not the other (e.g. a mis-declared plex or the wrong SDRF)

### Changed

- **BREAKING — Feature/PSM primary keys** (format 1.1, issue #217): feature PK
  `[sequence, charge, run_file_name, anchor_protein]` → `[peptidoform, charge,
  run_file_name, rt]`; psm PK `[sequence, charge, run_file_name, scan]` →
  `[peptidoform, charge, run_file_name, scan]`. Measured across ~13M real rows:
  `anchor_protein` is functionally redundant and apex `rt` is the only populated
  RT column that resolves co-eluting peaks of one peptidoform+charge in a run.
  Regenerate feature/psm files. `rt` is finite when present but **nullable**
  (some producers, e.g. FragPipe `combined_ion`, report no per-feature RT); a null
  `rt` is still a key value, so such producers must emit at most one feature per
  `(peptidoform, charge, run_file_name)` and any collision is caught as a
  `duplicate_pk` error; the key is within-file only.
- **Parquet output size**: writers now apply `BYTE_STREAM_SPLIT` encoding to high-entropy float columns (rt, rt_start/stop, predicted_rt, calculated/observed m/z, intensity arrays) and raise the ZSTD level to 9. Encoding-only and fully lossless — no schema change; output reads unchanged with pyarrow and DuckDB. Measured ~16% smaller on a 14 GB feature.parquet.

### Fixed

- **openms-consensus isobaric detection**: real quantms `IsobaricWorkflow` output stamps TMT/iTRAQ consensusXML with `experiment_type="label-free"` while the maps still carry `tmt6plex_*`/`itraq*plex_*` labels. The converter now detects channels from the **map label** (not `experiment_type`), so real quantms TMT no longer collapses all reporter channels into a single `LFQ` label. Verified on real cluster output (PXD000001 TMT → TMT126–131; BSA/PXD002395 LFQ unchanged)
- **RT unit conversion**: DIA-NN and MaxQuant converters now correctly convert retention time from minutes to seconds in feature and PSM parquet output
- **Code quality**: Spectronaut converter refactored to reduce cyclomatic complexity, fix logging f-string interpolation, remove unused arguments, and eliminate duplicate code
- **CDAP label-free intensity label**: label-free `PrecursorArea` intensities are now emitted with the `"LFQ"` label (aligned with the FragPipe/MaxQuant converters) so downstream label-free consumers (mokume) recognize them as primary intensities
- **QPX TMT/iTRAQ channel labels**: QuantMS mzTab feature `intensities[].label` now uses plex-aware canonical reporter names (TMT10 ch10 = `TMT131` not `TMT131N`, iTRAQ properly mapped); OpenMS `-out_qpx` enrichment relabels feature/pg intensities from filenames/bare indices to canonical names; `run.samples[].label` normalized to match — all from the shared sdrf-pipelines `channel_map` vocabulary

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
