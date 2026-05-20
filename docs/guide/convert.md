# Convert Commands

Convert various mass spectrometry data formats to the QPX standard format.

```python exec="1" session="doc_utils" result="ansi"
import click
import textwrap

def get_click_type_display(param):
    param_type = param.type
    type_str = str(param_type)
    if 'Path' in type_str:
        if hasattr(param_type, 'dir_okay') and not param_type.dir_okay:
            return 'FILE'
        elif hasattr(param_type, 'file_okay') and not param_type.file_okay:
            return 'DIRECTORY'
        else:
            return 'PATH'
    elif isinstance(param_type, click.types.FloatParamType):
        return 'FLOAT'
    elif isinstance(param_type, click.types.IntParamType):
        return 'INTEGER'
    elif param.is_flag:
        return 'FLAG'
    else:
        return 'TEXT'

def generate_params_table(command):
    table = '<table>\n<thead>\n<tr>\n'
    table += '<th>Parameter</th><th>Type</th><th>Required</th><th>Default</th><th>Description</th>\n'
    table += '</tr>\n</thead>\n<tbody>\n'
    for param in command.params:
        if isinstance(param, click.Option) and param.name not in ['help']:
            param_names = param.opts
            param_name = param_names[0] if param_names else f"--{param.name}"
            param_type = get_click_type_display(param)
            required = 'Yes' if param.required else 'No'

            # Extract default value from help text if it contains "(default: ...)"
            description = param.help or ''
            default_from_help = None
            if '(default:' in description.lower():
                import re
                match = re.search(r'\(default:\s*([^)]+)\)', description, re.IGNORECASE)
                if match:
                    default_from_help = match.group(1).strip()

            # Determine default value display
            if default_from_help:
                default = default_from_help
            elif param.default is not None and str(param.default) != 'Sentinel.UNSET':
                if param.is_flag:
                    default = '-'
                elif isinstance(param.default, (int, float)):
                    default = str(param.default)
                elif isinstance(param.default, str):
                    default = f'<code>{param.default}</code>'
                else:
                    default = str(param.default)
            else:
                default = '-'

            table += f'<tr>\n<td><code>{param_name}</code></td>\n<td>{param_type}</td>\n<td>{required}</td>\n<td>{default}</td>\n<td>{description}</td>\n</tr>\n'
    table += '</tbody>\n</table>'
    return table

def generate_description(command):
    if command.help:
        help_text = command.help
        if 'Example' in help_text:
            description = help_text.split('Example')[0].strip()
        else:
            description = help_text.strip()
        lines = description.split('\n')
        if len(lines) > 1:
            description = '\n'.join(lines[1:]).strip()
            return f'<p>{description}</p>'
    return ''

def generate_example(command, default_text=''):
    if command.help and 'Example' in command.help:
        example_section = command.help.split('Example')[1]
        if ':' in example_section:
            example_section = example_section.split(':', 1)[1]
        example_section = textwrap.dedent(example_section).strip()
        output = ''
        if default_text:
            output += f'<p>{default_text}</p>\n'
        output += f'<pre><code class="language-bash">{example_section}</code></pre>'
        return output
    return ''
```

## Overview

The `convert` command group provides converters for multiple proteomics software outputs, enabling standardization of data formats for downstream analysis. All commands generate parquet-format output files following the QPX specification.

## Available Commands

- [quantms](#quantms) - Convert QuantMS mzTab output to QPX format
- [diann](#diann) - Convert DIA-NN report to QPX format
- [spectronaut](#spectronaut) - Convert Spectronaut report to QPX format
- [maxquant](#maxquant) - Convert MaxQuant output to QPX format
- [fragpipe](#fragpipe) - Convert FragPipe output to QPX format
- [mzidentml](#mzidentml) - Convert mzIdentML file to PSM format
- [sdrf](#sdrf) - Convert SDRF to sample and run parquet files

---

## quantms

Convert QuantMS mzTab output to QPX format.

### Description {#quantms-description}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_quantms_cmd
print(generate_description(convert_quantms_cmd))
```

### Parameters {#quantms-parameters}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_quantms_cmd
print(generate_params_table(convert_quantms_cmd))
```

### Usage Examples {#quantms-examples}

#### Basic Example {#quantms-example-basic}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_quantms_cmd
print(generate_example(convert_quantms_cmd, 'Convert QuantMS data with default settings:'))
```

#### PSM Data Only {#quantms-example-psm}

```bash
qpxc convert quantms \
    --mztab-path tests/examples/quantms/dda-lfq-small/PXD007683-LFQ.sdrf_openms_design_openms.mzTab \
    --sdrf-file tests/examples/quantms/dda-lfq-small/PXD007683-LFQ.sdrf.tsv \
    --output-folder ./output \
    --structures psm \
    --output-prefix quantms_psm
```

#### Feature Data with MSstats {#quantms-example-feature}

```bash
qpxc convert quantms \
    --mztab-path tests/examples/quantms/dda-lfq-full/PXD007683-LFQ.sdrf_openms_design_openms.mzTab.gz \
    --msstats-file tests/examples/quantms/dda-lfq-full/PXD007683-LFQ.sdrf_openms_design_msstats_in.csv.gz \
    --sdrf-file tests/examples/quantms/dda-lfq-full/PXD007683-LFQ.sdrf.tsv \
    --output-folder ./output \
    --structures feature \
    --verbose
```

#### All Structures {#quantms-example-all}

```bash
qpxc convert quantms \
    --mztab-path tests/examples/quantms/dda-lfq-full/PXD007683-LFQ.sdrf_openms_design_openms.mzTab.gz \
    --msstats-file tests/examples/quantms/dda-lfq-full/PXD007683-LFQ.sdrf_openms_design_msstats_in.csv.gz \
    --sdrf-file tests/examples/quantms/dda-lfq-full/PXD007683-LFQ.sdrf.tsv \
    --output-folder ./output \
    --structures psm,feature,pg \
    --verbose
```

#### TMT Data {#quantms-example-tmt}

```bash
qpxc convert quantms \
    --mztab-path tests/examples/quantms/dda-plex-full/PXD007683TMT.sdrf_openms_design_openms.mzTab.gz \
    --msstats-file tests/examples/quantms/dda-plex-full/PXD007683TMT.sdrf_openms_design_msstats_in.csv.gz \
    --sdrf-file tests/examples/quantms/dda-plex-full/PXD007683-TMT.sdrf.tsv \
    --output-folder ./output \
    --structures pg
```

### Output Files {#quantms-output}

Depending on `--structures` parameter:

- **PSM**: `{output-prefix}-{uuid}.psm.parquet`
- **Feature**: `{output-prefix}-{uuid}.feature.parquet`
- **Protein Group**: `{output-prefix}-{uuid}.pg.parquet`

All files are in Parquet format and conform to their respective QPX specifications.

### Best Practices {#quantms-best-practices}

- Use `--structures` to control which output files are generated
- Provide `--msstats-file` when converting feature or pg structures
- Reuse database files with `--database-path` when processing the same mzTab multiple times
- Enable verbose mode for large datasets to monitor progress

---

## diann

Convert DIA-NN report files to QPX format.

### Description {#diann-description}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_diann_cmd
print(generate_description(convert_diann_cmd))
```

### Parameters {#diann-parameters}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_diann_cmd
print(generate_params_table(convert_diann_cmd))
```

### Usage Examples {#diann-examples}

#### Basic Example - Feature Data {#diann-example-basic}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_diann_cmd
print(generate_example(convert_diann_cmd, 'Convert a DIA-NN report with default settings:'))
```

#### Advanced Example with Partitioning {#diann-example-partitioning}

Convert with file partitioning based on run_file_name:

```bash
qpxc convert diann \
    --report-path tests/examples/diann/full/diann_report.tsv.gz \
    --qvalue-threshold 0.01 \
    --mzml-info-folder tests/examples/diann/full/mzml \
    --sdrf-file tests/examples/diann/full/PXD036609.sdrf.tsv \
    --output-folder ./output \
    --partitions run_file_name \
    --duckdb-max-memory 8GB \
    --duckdb-threads 4 \
    --verbose
```

#### Protein Groups from PG Matrix {#diann-example-pg}

Convert DIA-NN protein groups using the pg_matrix file:

```bash
qpxc convert diann \
    --report-path tests/examples/diann/full/diann_report.tsv.gz \
    --pg-matrix-path tests/examples/diann/full/diann_report.pg_matrix.tsv \
    --sdrf-file tests/examples/diann/full/PXD036609.sdrf.tsv \
    --output-folder ./output \
    --duckdb-max-memory 16GB \
    --duckdb-threads 8 \
    --verbose
```

### Output Files {#diann-output}

Output files generated:

- **Feature**: `{output-prefix}-{uuid}.feature.parquet` (always produced)
- **Protein Group**: `{output-prefix}-{uuid}.pg.parquet` (produced when `--pg-matrix-path` is provided)

### Common Issues {#diann-issues}

**Issue**: Out of memory errors with large files

- **Solution**: Increase `--duckdb-max-memory` parameter (e.g., `8GB`, `16GB`)

**Issue**: Slow processing

- **Solution**: Increase `--duckdb-threads` to utilize more CPU cores

**Issue**: Missing mzML info files

- **Solution**: Ensure all mzML info TSV files are in the specified folder with correct naming

### Best Practices {#diann-best-practices}

- Use Q-value threshold of 0.05 or lower for high-confidence results
- Enable partitioning for large datasets to improve memory usage
- Use verbose mode during initial testing to diagnose issues
- Ensure SDRF file correctly matches sample names in DIA-NN report
- For protein groups, ensure both report and pg_matrix files are from the same DIA-NN run

---

## spectronaut

Convert Spectronaut report files to QPX format.

### Description {#spectronaut-description}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_spectronaut_cmd
print(generate_description(convert_spectronaut_cmd))
```

### Parameters {#spectronaut-parameters}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_spectronaut_cmd
print(generate_params_table(convert_spectronaut_cmd))
```

### Usage Examples {#spectronaut-examples}

#### Basic Example {#spectronaut-example-basic}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_spectronaut_cmd
print(generate_example(convert_spectronaut_cmd, 'Convert Spectronaut data with default settings:'))
```

#### With SDRF Metadata {#spectronaut-example-sdrf}

```bash
qpxc convert spectronaut \
    --report-path report.tsv \
    --sdrf-file data.sdrf.tsv \
    --output-folder ./qpx_output \
    --qvalue-threshold 0.01 \
    --verbose
```

#### Large Dataset with Custom DuckDB Settings {#spectronaut-example-large}

```bash
qpxc convert spectronaut \
    --report-path large_report.tsv \
    --sdrf-file data.sdrf.tsv \
    --output-folder ./qpx_output \
    --duckdb-max-memory 16GB \
    --duckdb-threads 8 \
    --verbose
```

### Output Files {#spectronaut-output}

Output files generated:

- **Feature**: `{output-prefix}.feature.parquet` (always produced)
- **Protein Group**: `{output-prefix}.pg.parquet` (always produced)
- **Ontology**: `{output-prefix}.ontology.parquet` (score CV terms)
- **Provenance**: `{output-prefix}.provenance.parquet` (conversion metadata)

### Supported Spectronaut Columns {#spectronaut-columns}

The converter auto-detects column names from the Spectronaut report, supporting both standard and custom column naming. Key mapped columns include:

| Spectronaut Column       | QPX Field                  | Notes                              |
| ------------------------ | -------------------------- | ---------------------------------- |
| `R.FileName`             | `run_file_name`            | File extensions auto-stripped      |
| `PEP.StrippedSequence`   | `sequence`                 | Bare amino acid sequence           |
| `EG.ModifiedSequence`    | `peptidoform`              | Converted to ProForma notation     |
| `FG.Charge`              | `charge`                   | Precursor charge state             |
| `FG.Quantity`            | `intensities[0].intensity` | Raw precursor intensity            |
| `EG.Qvalue`              | `additional_scores`        | Mapped to PSI-MS CV terms          |
| `PG.ProteinGroups`       | `pg_accessions`            | Semicolon-delimited protein groups |
| `PG.Quantity`            | PG `intensities`           | Protein group quantity             |

### Common Issues {#spectronaut-issues}

**Issue**: European decimal separator (comma instead of period)

- **Solution**: The converter auto-detects decimal separators by inspecting the `EG.Qvalue` column

**Issue**: Out of memory with large reports

- **Solution**: Increase `--duckdb-max-memory` (e.g., `16GB`) and `--duckdb-threads`

**Issue**: Missing modifications in output

- **Solution**: Ensure `EG.ModifiedSequence` column is present in the Spectronaut export

### Best Practices {#spectronaut-best-practices}

- Export all required columns from Spectronaut (at minimum: `R.FileName`, `PEP.StrippedSequence`, `EG.ModifiedSequence`, `FG.Charge`, `FG.Quantity`, `PG.ProteinGroups`)
- Provide SDRF file for complete sample and run metadata
- Use Q-value threshold of 0.01 or 0.05 for high-confidence results
- For large datasets (>10M rows), increase DuckDB memory and threads

---

## maxquant

Convert MaxQuant output to QPX format.

### Description {#maxquant-description}

```python exec="1" session="doc_utils" html="1"
from qpx.cli.convert import convert_maxquant_cmd
print(generate_description(convert_maxquant_cmd))
```

### Parameters {#maxquant-parameters}

```python exec="1" session="doc_utils" html="1"
from qpx.cli.convert import convert_maxquant_cmd
print(generate_params_table(convert_maxquant_cmd))
```

### Usage Examples {#maxquant-examples}

#### Basic Example {#maxquant-example-basic}

```python exec="1" session="doc_utils" html="1"
from qpx.cli.convert import convert_maxquant_cmd
print(generate_example(convert_maxquant_cmd, 'Convert MaxQuant data with default settings:'))
```

#### PSM Data Only {#maxquant-example-psm}

```bash
qpxc convert maxquant \
    --msms-file tests/examples/maxquant/maxquant_simple/msms.txt \
    --output-folder ./output \
    --structures psm \
    --spectral-data \
    --output-prefix maxquant_psm
```

#### Feature Data with Protein Groups {#maxquant-example-feature}

```bash
qpxc convert maxquant \
    --evidence-file tests/examples/maxquant/maxquant_full/evidence.txt.gz \
    --protein-groups-file tests/examples/maxquant/maxquant_full/proteinGroups.txt \
    --sdrf-file tests/examples/maxquant/maxquant_full/PXD001819.sdrf.tsv \
    --output-folder ./output \
    --structures feature \
    --batch-size 500000 \
    --verbose
```

#### All Structures {#maxquant-example-all}

```bash
qpxc convert maxquant \
    --msms-file tests/examples/maxquant/maxquant_full/msms.txt.gz \
    --evidence-file tests/examples/maxquant/maxquant_full/evidence.txt.gz \
    --protein-groups-file tests/examples/maxquant/maxquant_full/proteinGroups.txt \
    --sdrf-file tests/examples/maxquant/maxquant_full/PXD001819.sdrf.tsv \
    --output-folder ./output \
    --structures psm,feature,pg \
    --batch-size 500000 \
    --verbose
```

### Output Files {#maxquant-output}

Depending on `--structures` parameter:

- **PSM**: `{output-prefix}-{uuid}.psm.parquet`
- **Feature**: `{output-prefix}-{uuid}.feature.parquet`
- **Protein Group**: `{output-prefix}-{uuid}.pg.parquet`

### Common Issues {#maxquant-issues}

**Issue**: Memory errors with compressed evidence files

- **Solution**: Reduce `--batch-size` or increase available RAM

**Issue**: Missing Q-value information

- **Solution**: Provide `--protein-groups-file` for accurate Q-value mapping

### Best Practices {#maxquant-best-practices}

- Use `--structures` to control which output files are generated
- Always provide `--protein-groups-file` when available for better data quality
- Ensure SDRF sample names match MaxQuant experiment names
- Use compressed files (.gz) to save disk space
- Adjust `--batch-size` based on available memory
- Use `--spectral-data` flag if downstream analysis requires spectral information

---

## fragpipe

Convert FragPipe output to QPX format.

### Description {#fragpipe-description}

```python exec="1" session="doc_utils" html="1"
from qpx.cli.convert import convert_fragpipe_cmd
print(generate_description(convert_fragpipe_cmd))
```

### Parameters {#fragpipe-parameters}

```python exec="1" session="doc_utils" html="1"
from qpx.cli.convert import convert_fragpipe_cmd
print(generate_params_table(convert_fragpipe_cmd))
```

### Usage Examples {#fragpipe-examples}

#### Basic Example {#fragpipe-example-basic}

```python exec="1" session="doc_utils" html="1"
from qpx.cli.convert import convert_fragpipe_cmd
print(generate_example(convert_fragpipe_cmd, 'Convert FragPipe PSM data with default settings:'))
```

#### With Custom Settings {#fragpipe-example-custom}

```bash
qpxc convert fragpipe \
    --msms-file /path/to/psm.tsv \
    --output-folder ./output \
    --batch-size 500000 \
    --output-prefix fragpipe_psm
```

### Output Files {#fragpipe-output}

- **Output**: `{output-prefix}-{uuid}.psm.parquet`
- **Format**: Parquet file containing PSM data
- **Schema**: Conforms to QPX PSM specification

---

## mzidentml

Convert mzIdentML (.mzid) files to QPX PSM parquet format.

### Description {#mzidentml-description}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_mzidentml_cmd
print(generate_description(convert_mzidentml_cmd))
```

### Parameters {#mzidentml-parameters}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_mzidentml_cmd
print(generate_params_table(convert_mzidentml_cmd))
```

### Usage Examples {#mzidentml-examples}

#### Basic Example {#mzidentml-example-basic}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_mzidentml_cmd
print(generate_example(convert_mzidentml_cmd, 'Convert an mzIdentML file with default settings:'))
```

#### With Spectral Data from Single mzML {#mzidentml-example-spectral}

```bash
qpxc convert mzidentml \
    --mzid-file /path/to/results.mzid \
    --mzml-file /path/to/spectra.mzML \
    --output-folder ./output \
    --spectral-data \
    --output-prefix psm_with_spectra
```

#### With Spectral Data from Multiple mzML Files {#mzidentml-example-multi-mzml}

When your mzIdentML references multiple mzML files, use the `--mzml-folder` option:

```bash
qpxc convert mzidentml \
    --mzid-file /path/to/results.mzid.gz \
    --mzml-folder /path/to/mzml_files/ \
    --output-folder ./output \
    --spectral-data \
    --output-prefix psm_multi_mzml
```

The converter automatically matches PSMs to mzML files based on the `run_file_name` field in the mzIdentML. File matching is case-insensitive and supports both `.mzML` and `.mzML.gz` extensions.

### Supported Native ID Formats {#mzidentml-native-id}

The converter supports multiple native ID formats for scan number extraction:

| Format | Vendor/Source | Example |
|--------|---------------|---------|
| `scan=XXX` | Thermo | `controllerType=0 controllerNumber=1 scan=12345` |
| `cycle=XXX` | Waters/Agilent | `sample=1 period=1 cycle=1055 experiment=4` |
| `index=XXX` | Generic | `index=500` |
| `spectrum=XXX` | Various | `spectrum=999` |

### Output Files {#mzidentml-output}

- **Output**: `{output-prefix}-{uuid}.psm.parquet`
- **Format**: Parquet file containing PSM-level data
- **Schema**: Conforms to QPX PSM specification

### Supported mzIdentML Features {#mzidentml-features}

- **Compressed files**: Supports both `.mzid` and `.mzid.gz` formats
- **Modifications**: Full support for UNIMOD and custom modifications
- **Scores**: Extracts all CV-term scores with `higher_better` flag annotation
- **Decoy detection**: Automatic detection via `isDecoy` attribute
- **Multi-file support**: Handles mzIdentML referencing multiple spectra files

### Best Practices {#mzidentml-best-practices}

- Use `--mzml-folder` when mzIdentML references multiple mzML files
- Ensure mzML file names match those referenced in mzIdentML (case-insensitive)
- Use compressed `.mzid.gz` files to save disk space
- Enable `--spectral-data` only when spectral arrays are needed for downstream analysis

### Common Issues {#mzidentml-issues}

**Issue**: No spectra attached from mzML folder

- **Solution**: Verify mzML file names match `run_file_name` in mzIdentML

**Issue**: zlib errors when reading mzML.gz files

- **Solution**: Decompress mzML.gz files or re-download if corrupted

**Issue**: Scan numbers not extracted correctly

- **Solution**: Check if your native ID format is supported; the converter auto-detects common formats

---

## sdrf

Convert SDRF metadata files to QPX sample and run parquet format.

### Description {#sdrf-description}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_sdrf_cmd
print(generate_description(convert_sdrf_cmd))
```

### Parameters {#sdrf-parameters}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_sdrf_cmd
print(generate_params_table(convert_sdrf_cmd))
```

### Usage Examples {#sdrf-examples}

#### Basic Example {#sdrf-example-basic}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.convert import convert_sdrf_cmd
print(generate_example(convert_sdrf_cmd, 'Convert SDRF metadata with default settings:'))
```

### Output Files {#sdrf-output}

- **Sample**: `{output-prefix}-{uuid}.sample.parquet`
- **Run**: `{output-prefix}-{uuid}.run.parquet`
- **Format**: Parquet files containing sample and run metadata
- **Schema**: Conforms to QPX sample and run specifications

### Best Practices {#sdrf-best-practices}

- Ensure SDRF file follows the PRIDE SDRF specifications
- Use verbose mode to diagnose parsing issues
- The converter automatically maps SDRF characteristics to QPX ontology terms

---

## Related Commands

- [Transform Commands](transform.md) - Further process converted data
- [Visualization Commands](visualize.md) - Create plots from converted data
- [Statistics Commands](stats.md) - Generate statistics from converted data
