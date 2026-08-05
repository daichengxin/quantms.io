# Converter Coverage Matrix

This page lists which QPX data views each converter produces. Use it to see at a glance what output to expect from each tool.

## Views produced per converter

| Converter | PSM | Feature | PG | Pepmap | Sample | Run | Dataset | Ontology | Provenance | mz |
|-----------|:---:|:-------:|:--:|:------:|:------:|:---:|:-------:|:--------:|:----------:|:--:|
| **MaxQuant** | Yes | Yes | Yes | No | If SDRF | If SDRF | Yes | If SDRF | No | No |
| **FragPipe** | Yes | Yes | Yes | No | If SDRF | If SDRF | Yes | If SDRF | No | No |
| **DIA-NN** | No | Yes | Yes | No | If SDRF | If SDRF | Yes | If SDRF | No | No |
| **Spectronaut** | No | Yes | Yes | No | If SDRF | If SDRF | Yes | Yes | Yes | No |
| **OpenMS native QPX** | Yes | Yes | Yes | No | Yes | Yes | Yes | Yes | Yes | No |
| **OpenMS consensusXML** | Yes | Yes | Yes | No | If SDRF | If SDRF | No | No | No | No |
| **CDAP** | Yes | Yes | Yes | No | If PDC | If PDC | Yes | Yes | Yes | No |
| **mzIdentML** | Yes | No | No | Yes | If SDRF | If SDRF | Yes | Yes | Yes | No |
| **SDRF** | No | No | No | No | Yes | Yes | No | Optional | No | No |

- **Yes** — the converter produces this view.
- **No** — the converter does not produce this view (e.g. DIA-NN has no PSM view; mzIdentML has no Feature/PG).
- **If SDRF** — the view is produced only when an SDRF file is provided (sample and run metadata). SDRF is optional except for the `qpxc convert openms` command, which requires it.
- **If PDC** — CPTAC/PDC studies ship no SDRF, and CDAP `.psm` files carry no sample metadata. When run through `qpxc pdc2qpx` (default), the sample/run views are built from PDC GraphQL metadata, which also recovers the TMT/iTRAQ channel → biological-sample mapping. Disable with `--no-metadata`.

> The `mz` (full-spectra) view is produced by the standalone `qpxc convert mz` command, or automatically by `qpxc pdc2qpx --include-spectra`; it is not emitted by the per-tool converters above.

## CLI commands

| Converter | Command |
|-----------|---------|
| MaxQuant | `qpxc convert maxquant` |
| FragPipe | `qpxc convert fragpipe` |
| DIA-NN | `qpxc convert diann` |
| Spectronaut | `qpxc convert spectronaut` |
| OpenMS native QPX | `qpxc convert openms` |
| OpenMS consensusXML | `qpxc convert openms-consensus` |
| CDAP | `qpxc convert cdap` |
| mzIdentML | `qpxc convert mzidentml` |
| SDRF only | `qpxc convert sdrf` |

## Input files (summary)

| Converter | Typical inputs |
|-----------|----------------|
| MaxQuant | msms.txt, evidence.txt, proteinGroups.txt |
| FragPipe | psm.tsv, combined_ion/combined_peptide, combined_protein |
| DIA-NN | report (tsv), pg_matrix (optional) |
| Spectronaut | report.tsv; optional SDRF |
| OpenMS native QPX | `-out_qpx` Parquet directory; SDRF; optional companion consensusXML |
| OpenMS consensusXML | consensusXML; optional SDRF |
| CDAP | CPTAC CDAP `.psm` files in one study directory |
| mzIdentML | .mzid / .mzid.gz; optional MGF or mzML (file/folder) for spectra; optional SDRF |
| SDRF | Single SDRF TSV file |

For field-level mappings from each tool’s columns to QPX, see [Tool Field Mappings](tool-mappings.md).
