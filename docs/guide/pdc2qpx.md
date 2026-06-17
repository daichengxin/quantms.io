# pdc2qpx

One-shot download + conversion of a PDC/CPTAC study into a QPX dataset.

`pdc2qpx` orchestrates three steps in a single command:

1. Download the study's CDAP `.psm` files (and, with `--include-spectra`, the mzML spectra) from PDC via [pridepy](https://github.com/bigbio/pridepy).
2. Convert the `.psm` files to QPX psm/feature/pg/dataset/ontology/provenance views (`CdapConverter`).
3. Convert the mzML files to a full-spectra `mz.parquet`, linked to PSM/feature by `run_file_name` + `scan`.

The download layer requires the optional `pdc` extra:

```bash
pip install qpx[pdc]
```

## Usage

```bash
qpxc pdc2qpx --study PDC000109 --download-dir ./downloads --output-folder ./qpx/PDC000109 [OPTIONS]
```

## Parameters

| Option | Description | Default |
| ------ | ----------- | ------- |
| `--study` | PDC study accession (e.g. `PDC000109`) | required |
| `--download-dir` | Directory for downloaded files (written under `<download-dir>/<study>/`) | required |
| `--output-folder` | Directory for generated QPX parquet files | required |
| `--include-spectra` | Also download mzML and produce a full-spectra `<study>.mz.parquet` | off |
| `--ms-levels` | MS levels for the mz view (e.g. `2` or `1,2`) | all |
| `--max-cpus` | Threads for the CDAP conversion | 24 |
| `--max-memory` | Memory limit for the CDAP conversion | 16GB |
| `--download-threads` | Parallel HTTP Range threads per file | 24 |
| `--skip-download` | Use files already present under `<download-dir>/<study>/` | off |

## Examples

### Base QPX (psm/feature/pg) only

```bash
qpxc pdc2qpx \
    --study PDC000109 \
    --download-dir ./downloads \
    --output-folder ./qpx/PDC000109
```

### Entire QPX including full spectra

For quantms-style reanalysis you need MS1 (precursor LFQ quantification) as well as MS2 (identification), so include all MS levels (the default when `--ms-levels` is omitted):

```bash
qpxc pdc2qpx \
    --study PDC000109 \
    --download-dir ./downloads \
    --output-folder ./qpx/PDC000109 \
    --include-spectra --max-cpus 24
```

### Re-run on already-downloaded files

```bash
qpxc pdc2qpx \
    --study PDC000109 \
    --download-dir ./downloads \
    --output-folder ./qpx/PDC000109 \
    --include-spectra --skip-download
```

## Output

Under `--output-folder`, prefixed with the study accession:

- `<study>.psm.parquet`, `<study>.feature.parquet`, `<study>.pg.parquet`
- `<study>.dataset.parquet`, `<study>.ontology.parquet`, `<study>.provenance.parquet`
- `<study>.mz.parquet` (only with `--include-spectra`) — full spectra, joinable to PSM/feature by `run_file_name` + `scan`

## Notes

- Only CPTAC studies that publish CDAP `.psm` files can produce the base QPX views; `--include-spectra` additionally requires mzML to be available on PDC for that study.
- Downloads go through pridepy's PDC layer (signed CloudFront URLs, md5 validation, automatic 403 URL re-signing).
- See [Converters](convert.md) for the underlying `qpxc convert cdap` and `qpxc convert mz` commands.
