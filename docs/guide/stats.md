# Statistics and Quality Control

Generate statistics and perform quality control analysis on QPX datasets using the Python API.

## Overview

QPX provides a comprehensive Python API for accessing dataset statistics and quality control metrics. The API allows you to programmatically query your data and generate statistical summaries for publications, quality assessment, and comparative analysis.

## Getting Started

Load a QPX dataset and access basic statistics:

```python
import qpx

# Load dataset
ds = qpx.Dataset("path/to/dataset/")

# Basic counts
print(f"Total PSMs: {ds.psm.count()}")
print(f"Total features: {ds.feature.count()}")
print(f"Total protein groups: {ds.pg.count()}")
```

## Core Statistics

### PSM-Level Statistics

Access peptide-spectrum match statistics:

```python
import qpx

ds = qpx.Dataset("./output/")

# PSM counts
total_psms = ds.psm.count()
print(f"Number of PSMs: {total_psms:,}")

# Unique peptides (without modifications)
unique_peptides = ds.psm.data['sequence'].nunique()
print(f"Number of peptides: {unique_peptides:,}")

# Unique peptidoforms (with modifications)
unique_peptidoforms = ds.psm.data['peptidoform'].nunique()
print(f"Number of peptidoforms: {unique_peptidoforms:,}")

# Unique proteins
unique_proteins = ds.psm.data['protein_accessions'].nunique()
print(f"Number of proteins: {unique_proteins:,}")

# Number of MS runs
num_runs = ds.psm.data['run_file_name'].nunique()
print(f"Number of msruns: {num_runs:,}")
```

### Feature-Level Statistics

Access feature quantification statistics:

```python
import qpx

ds = qpx.Dataset("./output/")

# Feature counts
total_features = ds.feature.count()
print(f"Total features: {total_features:,}")

# Samples with quantification data
num_samples = len([col for col in ds.feature.data.columns
                   if col.startswith('sample_')])
print(f"Number of samples: {num_samples}")

# Missing value analysis
intensity_cols = [col for col in ds.feature.data.columns
                 if col.startswith('sample_')]
missing_values = ds.feature.data[intensity_cols].isna().sum().sum()
total_values = ds.feature.data[intensity_cols].size
completeness = (1 - missing_values / total_values) * 100
print(f"Data completeness: {completeness:.2f}%")
```

### Protein Group Statistics

Access protein group quantification statistics:

```python
import qpx

ds = qpx.Dataset("./output/")

# Protein group counts
total_pgs = ds.pg.count()
print(f"Number of protein groups: {total_pgs:,}")

# iBAQ statistics (if available)
if 'ibaq' in ds.pg.data.columns:
    ibaq_cols = [col for col in ds.pg.data.columns
                 if col.startswith('sample_') and 'ibaq' in col.lower()]
    num_ibaq_samples = len(ibaq_cols)
    num_ibaq_proteins = ds.pg.data[ibaq_cols].notna().any(axis=1).sum()
    print(f"iBAQ Number of proteins: {num_ibaq_proteins:,}")
    print(f"iBAQ Number of samples: {num_ibaq_samples}")
```

## Comprehensive Statistics Report

Generate a complete statistics summary:

```python
import qpx

def generate_stats_report(dataset_path, output_file=None):
    """Generate comprehensive statistics for a QPX dataset."""
    ds = qpx.Dataset(dataset_path)

    report = []
    report.append("=" * 50)
    report.append("QPX Dataset Statistics Report")
    report.append("=" * 50)
    report.append("")

    # PSM Statistics
    if hasattr(ds, 'psm') and ds.psm.count() > 0:
        report.append("PSM Statistics:")
        report.append(f"  Number of PSMs: {ds.psm.count():,}")
        report.append(f"  Number of peptides: {ds.psm.data['sequence'].nunique():,}")
        report.append(f"  Number of peptidoforms: {ds.psm.data['peptidoform'].nunique():,}")
        report.append(f"  Number of proteins: {ds.psm.data['protein_accessions'].nunique():,}")
        report.append(f"  Number of msruns: {ds.psm.data['run_file_name'].nunique():,}")
        report.append("")

    # Feature Statistics
    if hasattr(ds, 'feature') and ds.feature.count() > 0:
        intensity_cols = [col for col in ds.feature.data.columns
                         if col.startswith('sample_')]
        report.append("Feature Statistics:")
        report.append(f"  Number of features: {ds.feature.count():,}")
        report.append(f"  Number of samples: {len(intensity_cols)}")
        report.append("")

    # Protein Group Statistics
    if hasattr(ds, 'pg') and ds.pg.count() > 0:
        report.append("Protein Group Statistics:")
        report.append(f"  Number of protein groups: {ds.pg.count():,}")
        report.append("")

    # Dataset metadata
    if hasattr(ds, 'dataset_meta'):
        report.append("Dataset Metadata:")
        if hasattr(ds.dataset_meta, 'name'):
            report.append(f"  Name: {ds.dataset_meta.name}")
        if hasattr(ds.dataset_meta, 'version'):
            report.append(f"  Version: {ds.dataset_meta.version}")

    report.append("=" * 50)

    # Output
    report_text = "\n".join(report)
    if output_file:
        with open(output_file, 'w') as f:
            f.write(report_text)
        print(f"Report saved to: {output_file}")
    else:
        print(report_text)

    return report_text

# Usage
generate_stats_report("./output/", "./reports/statistics.txt")
```

## Statistical Interpretation Guide

### Data Quality Indicators

#### Good Quality Signs

- Consistent protein counts across replicates
- Expected PSM-to-peptide ratios (typically 2-5)
- Complete data across all msruns
- Reasonable peptidoform diversity (typically 1-3 per peptide)

#### Potential Issues

- **Very low PSM counts**: Search parameter issues, poor sample quality
- **High peptidoform-to-peptide ratio**: Over-prediction of modifications
- **Missing msruns**: File processing errors
- **Extreme variation**: Batch effects, contamination

### Understanding Metrics

#### Peptide vs Peptidoform

- **Peptide**: Amino acid sequence (e.g., `PEPTIDE`)
- **Peptidoform**: Peptide + modifications (e.g., `PEPTIDE[+16]`)

#### Expected Ratios

Typical ratios for quality data:

- **PSMs per peptide**: 2-5 (varies with replicates)
- **Peptidoforms per peptide**: 1-3 (depends on PTM analysis)
- **Peptides per protein**: 3-20 (depends on protein abundance and coverage)

### Comparative Analysis

When comparing datasets:

1. **Normalize by sample amount**: Account for loading differences
2. **Consider instrument type**: Different platforms yield different numbers
3. **Account for search space**: Database size affects identification rates
4. **Match FDR thresholds**: Ensure fair comparison

### Reporting Guidelines

For publications, report:

- Total unique proteins, peptides, and PSMs
- Number of biological replicates
- Number of technical replicates (msruns)
- FDR thresholds applied
- Protein and peptide identification rates

## Quality Control Workflow

Complete QC workflow example:

```python
import qpx
import pandas as pd

def quality_control_report(dataset_path):
    """Perform comprehensive quality control analysis."""
    ds = qpx.Dataset(dataset_path)

    # 1. Check data completeness
    print("1. Data Completeness Check")
    print(f"   PSM data: {'Present' if ds.psm.count() > 0 else 'Missing'}")
    print(f"   Feature data: {'Present' if ds.feature.count() > 0 else 'Missing'}")
    print(f"   Protein group data: {'Present' if ds.pg.count() > 0 else 'Missing'}")
    print()

    # 2. Identification rates
    print("2. Identification Rates")
    total_psms = ds.psm.count()
    unique_peptides = ds.psm.data['sequence'].nunique()
    unique_proteins = ds.psm.data['protein_accessions'].nunique()
    print(f"   PSMs per peptide: {total_psms / unique_peptides:.2f}")
    print(f"   Peptides per protein: {unique_peptides / unique_proteins:.2f}")
    print()

    # 3. Missing value analysis
    if hasattr(ds, 'feature') and ds.feature.count() > 0:
        print("3. Missing Value Analysis")
        intensity_cols = [col for col in ds.feature.data.columns
                         if col.startswith('sample_')]
        missing_pct = ds.feature.data[intensity_cols].isna().mean() * 100
        print(f"   Average missing values: {missing_pct.mean():.2f}%")
        print(f"   Range: {missing_pct.min():.2f}% - {missing_pct.max():.2f}%")
        print()

    # 4. Modification analysis
    print("4. Modification Analysis")
    peptidoforms = ds.psm.data['peptidoform'].nunique()
    peptides = ds.psm.data['sequence'].nunique()
    mod_ratio = peptidoforms / peptides
    print(f"   Peptidoforms per peptide: {mod_ratio:.2f}")
    print(f"   Status: {'Normal' if 1 <= mod_ratio <= 3 else 'Review recommended'}")
    print()

    # 5. Run completeness
    print("5. MS Run Completeness")
    num_runs = ds.psm.data['run_file_name'].nunique()
    print(f"   Total runs: {num_runs}")
    run_psms = ds.psm.data.groupby('run_file_name').size()
    print(f"   PSMs per run: {run_psms.mean():.0f} ± {run_psms.std():.0f}")
    print()

    return {
        'total_psms': total_psms,
        'unique_peptides': unique_peptides,
        'unique_proteins': unique_proteins,
        'num_runs': num_runs,
        'mod_ratio': mod_ratio
    }

# Usage
qc_results = quality_control_report("./output/")
```

## Custom Thresholds

Define and check quality thresholds:

```python
import qpx

def check_quality_thresholds(dataset_path):
    """Check dataset against defined quality thresholds."""
    ds = qpx.Dataset(dataset_path)

    # Define thresholds
    thresholds = {
        'min_proteins': 1000,
        'min_peptides': 5000,
        'min_psms': 10000,
        'max_missing_pct': 30,
        'min_runs': 3
    }

    # Check against thresholds
    warnings = []

    # Protein count
    protein_count = ds.psm.data['protein_accessions'].nunique()
    if protein_count < thresholds['min_proteins']:
        warnings.append(f"Low protein count: {protein_count:,} (threshold: {thresholds['min_proteins']:,})")

    # Peptide count
    peptide_count = ds.psm.data['sequence'].nunique()
    if peptide_count < thresholds['min_peptides']:
        warnings.append(f"Low peptide count: {peptide_count:,} (threshold: {thresholds['min_peptides']:,})")

    # PSM count
    psm_count = ds.psm.count()
    if psm_count < thresholds['min_psms']:
        warnings.append(f"Low PSM count: {psm_count:,} (threshold: {thresholds['min_psms']:,})")

    # Run count
    run_count = ds.psm.data['run_file_name'].nunique()
    if run_count < thresholds['min_runs']:
        warnings.append(f"Low run count: {run_count} (threshold: {thresholds['min_runs']})")

    # Report
    if warnings:
        print("Quality Warnings:")
        for warning in warnings:
            print(f"  - {warning}")
    else:
        print("All quality thresholds passed!")

    return len(warnings) == 0

# Usage
passes_qc = check_quality_thresholds("./output/")
```

## Use Cases

### Quality Control

Verify expected number of identifications:

```python
import qpx

ds = qpx.Dataset("./output/")
stats = {
    'proteins': ds.psm.data['protein_accessions'].nunique(),
    'peptides': ds.psm.data['sequence'].nunique(),
    'psms': ds.psm.count()
}

print(f"Identified {stats['proteins']:,} proteins from {stats['psms']:,} PSMs")
```

### Publication Reporting

Generate summary statistics for methods sections:

```python
import qpx

ds = qpx.Dataset("./output/")

print("Statistical Summary for Publication:")
print(f"  - {ds.psm.data['protein_accessions'].nunique():,} unique proteins")
print(f"  - {ds.psm.data['sequence'].nunique():,} unique peptides")
print(f"  - {ds.psm.count():,} peptide-spectrum matches")
print(f"  - {ds.psm.data['run_file_name'].nunique()} MS runs")
```

### Data Completeness

Assess coverage across samples:

```python
import qpx

ds = qpx.Dataset("./output/")

if hasattr(ds, 'feature') and ds.feature.count() > 0:
    intensity_cols = [col for col in ds.feature.data.columns
                     if col.startswith('sample_')]

    # Per-sample completeness
    for col in intensity_cols:
        completeness = (1 - ds.feature.data[col].isna().mean()) * 100
        print(f"{col}: {completeness:.1f}% complete")
```

### Method Optimization

Compare different search parameters:

```python
import qpx

# Load two datasets processed with different parameters
ds1 = qpx.Dataset("./output_params1/")
ds2 = qpx.Dataset("./output_params2/")

print("Parameter Set 1:")
print(f"  Proteins: {ds1.psm.data['protein_accessions'].nunique():,}")
print(f"  Peptides: {ds1.psm.data['sequence'].nunique():,}")
print()

print("Parameter Set 2:")
print(f"  Proteins: {ds2.psm.data['protein_accessions'].nunique():,}")
print(f"  Peptides: {ds2.psm.data['sequence'].nunique():,}")
```

## Batch Processing

Process multiple datasets:

```python
import qpx
from pathlib import Path

def batch_statistics(dataset_paths, output_file):
    """Generate statistics for multiple datasets."""
    results = []

    for path in dataset_paths:
        ds = qpx.Dataset(path)
        stats = {
            'dataset': Path(path).name,
            'proteins': ds.psm.data['protein_accessions'].nunique(),
            'peptides': ds.psm.data['sequence'].nunique(),
            'psms': ds.psm.count(),
            'runs': ds.psm.data['run_file_name'].nunique()
        }
        results.append(stats)

    # Save to file
    import pandas as pd
    df = pd.DataFrame(results)
    df.to_csv(output_file, index=False)
    print(f"Batch statistics saved to: {output_file}")

    return df

# Usage
datasets = ["./dataset1/", "./dataset2/", "./dataset3/"]
batch_statistics(datasets, "./reports/batch_stats.csv")
```

## Related Documentation

- [Visualization Guide](visualize.md) - Create visual representations of statistics
- [Dataset API Reference](../spec/dataset.md) - Core dataset structure
- [PSM Format](../spec/psm.md) - PSM data specification
- [Feature Format](../spec/feature.md) - Feature data specification
- [Protein Group Format](../spec/pg.md) - Protein group data specification
