# Visualization Guide

Create publication-quality visualizations from QPX datasets using the Python API.

## Overview

QPX provides a comprehensive visualization API through dataset views. Each view has a `.plot()` method that generates informative, publication-ready plots. All plots are created using matplotlib and can be saved in vector formats (SVG, PDF) for high-resolution output.

## Getting Started

Load a QPX dataset and create visualizations:

```python
import qpx

# Load dataset
ds = qpx.Dataset("path/to/dataset/")

# View-based plotting
ds.identifications.plot()  # Identification summary plot
ds.runs.plot()             # Run-level QC plot
ds.modifications.plot()    # Modification distribution plot
ds.qc.plot()              # Quality control dashboard
```

## Available Visualizations

### Identification Summary

Visualize protein, peptide, and PSM identifications:

```python
import qpx

ds = qpx.Dataset("./output/")

# Generate identification summary plot
fig = ds.identifications.plot()

# Save to file
fig.savefig("./plots/identifications.svg", format='svg', bbox_inches='tight')
```

**Output**:

- Bar plot showing counts of proteins, peptides, and PSMs
- Stacked bars for different identification levels
- Useful for quick assessment of dataset size

**Interpretation**:

- **High protein counts**: Good depth of coverage
- **Low peptide-to-protein ratio**: May indicate poor fragmentation or search issues
- **High PSM-to-peptide ratio**: Good reproducibility across runs

### Run Summary

Visualize statistics across MS runs:

```python
import qpx

ds = qpx.Dataset("./output/")

# Generate run-level QC plot
fig = ds.runs.plot()

# Save to file
fig.savefig("./plots/run_summary.svg", format='svg', bbox_inches='tight')
```

**Output**:

- Multiple panels showing run-level metrics
- PSM counts per run
- Peptide and protein identifications per run
- Useful for identifying problematic runs

**Interpretation**:

- **Consistent bars**: Good technical reproducibility
- **Outlier runs**: May indicate instrument issues or sample problems
- **Declining counts**: Possible column degradation

### Modification Distribution

Visualize post-translational modifications:

```python
import qpx

ds = qpx.Dataset("./output/")

# Generate modification distribution plot
fig = ds.modifications.plot()

# Save to file
fig.savefig("./plots/modifications.svg", format='svg', bbox_inches='tight')
```

**Output**:

- Bar plot showing frequency of different modifications
- Grouped by modification type
- Useful for PTM analysis validation

**Interpretation**:

- **Expected modifications**: Oxidation, carbamidomethylation, etc.
- **Unexpected modifications**: May indicate search parameter issues
- **Modification frequency**: Reflects biological state and search sensitivity

### Quality Control Dashboard

Comprehensive QC visualization:

```python
import qpx

ds = qpx.Dataset("./output/")

# Generate QC dashboard
fig = ds.qc.plot()

# Save to file
fig.savefig("./plots/qc_dashboard.svg", format='svg', bbox_inches='tight')
```

**Output**:

- Multi-panel dashboard with key QC metrics
- Intensity distributions
- Missing value patterns
- Identification rates
- Run-to-run consistency

**Interpretation**:

- **Aligned intensity distributions**: Good normalization
- **Low missing values**: Complete quantification
- **Consistent identification rates**: Technical reproducibility

## Intensity Distribution Plots

### Box Plot

Visualize intensity distributions across samples:

```python
import qpx
import matplotlib.pyplot as plt

ds = qpx.Dataset("./output/")

# Get intensity columns
intensity_cols = [col for col in ds.feature.data.columns
                 if col.startswith('sample_')]

# Create box plot
fig, ax = plt.subplots(figsize=(12, 6))
ds.feature.data[intensity_cols].apply(lambda x: x[x > 0].apply('log10')).boxplot(ax=ax)
ax.set_xlabel('Sample')
ax.set_ylabel('log10(Intensity)')
ax.set_title('Intensity Distribution Across Samples')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()

# Save
fig.savefig("./plots/intensity_boxplot.svg", format='svg', bbox_inches='tight')
```

**Output**:

- Box plots for each sample showing intensity distribution
- Log-transformed intensities for better visualization
- Box shows interquartile range (IQR)
- Whiskers extend to 1.5×IQR
- Outliers shown as individual points

**Interpretation**:

- **Aligned medians**: Good normalization
- **Similar IQR**: Consistent quantification across samples
- **Many outliers**: May indicate contamination or technical issues
- **Different ranges**: Batch effects or loading differences

### KDE (Kernel Density Estimation)

Plot smooth density distributions:

```python
import qpx
import matplotlib.pyplot as plt
import numpy as np

ds = qpx.Dataset("./output/")

# Get intensity columns (limit to first 10 samples for readability)
intensity_cols = [col for col in ds.feature.data.columns
                 if col.startswith('sample_')][:10]

# Create KDE plot
fig, ax = plt.subplots(figsize=(10, 6))
for col in intensity_cols:
    intensities = ds.feature.data[col].dropna()
    if len(intensities) > 0:
        log_intensities = np.log10(intensities[intensities > 0])
        log_intensities.plot.kde(ax=ax, label=col)

ax.set_xlabel('log10(Intensity)')
ax.set_ylabel('Density')
ax.set_title('Intensity Distribution (KDE)')
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save
fig.savefig("./plots/intensity_kde.svg", format='svg', bbox_inches='tight')
```

**Output**:

- Overlaid kernel density curves for each sample
- Smooth representation of intensity distributions
- Legend shows sample identifiers

**Interpretation**:

- **Overlapping curves**: Good sample-to-sample consistency
- **Shifted curves**: Potential batch effects or normalization issues
- **Different shapes**: Sample-specific technical issues
- **Bimodal distributions**: Distinct protein abundance classes

## Peptide Distribution Plots

### Peptides per Protein

Visualize peptide coverage across proteins:

```python
import qpx
import matplotlib.pyplot as plt

ds = qpx.Dataset("./output/")

# Count peptides per protein
peptides_per_protein = ds.psm.data.groupby('protein_accessions')['sequence'].nunique()
top_proteins = peptides_per_protein.nlargest(20)

# Create bar plot
fig, ax = plt.subplots(figsize=(12, 6))
top_proteins.plot.bar(ax=ax)
ax.set_xlabel('Protein')
ax.set_ylabel('Number of Peptides')
ax.set_title('Peptide Distribution Across Top 20 Proteins')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()

# Save
fig.savefig("./plots/peptides_per_protein.svg", format='svg', bbox_inches='tight')
```

**Output**:

- Bar plot showing peptide counts per protein
- Top N proteins by peptide count
- X-axis: Protein identifiers
- Y-axis: Number of unique peptides

**Interpretation**:

- **High peptide counts**: Abundant proteins with good coverage
- **Single peptide proteins**: May be less confident identifications
- **Distribution shape**: Reflects proteome complexity

### Peptides by Condition

Compare peptide identifications across experimental conditions:

```python
import qpx
import matplotlib.pyplot as plt

ds = qpx.Dataset("./output/")

# Assuming 'condition' column exists in PSM data or can be derived from sample metadata
# This is a simplified example - adapt based on your metadata structure
if 'condition' in ds.psm.data.columns:
    peptides_by_condition = ds.psm.data.groupby('condition')['sequence'].nunique()

    # Create bar plot
    fig, ax = plt.subplots(figsize=(10, 6))
    peptides_by_condition.plot.bar(ax=ax)
    ax.set_xlabel('Condition')
    ax.set_ylabel('Number of Peptides')
    ax.set_title('Peptides Identified by Condition')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    # Save
    fig.savefig("./plots/peptides_by_condition.svg", format='svg', bbox_inches='tight')
```

**Interpretation**:

- **High variation**: May indicate batch effects or quality issues
- **Low counts in specific conditions**: May suggest technical problems
- **Consistent counts**: Good data quality and reproducibility

## iBAQ Distribution

Visualize absolute protein abundance (iBAQ values):

```python
import qpx
import matplotlib.pyplot as plt
import numpy as np

ds = qpx.Dataset("./output/")

# Assuming iBAQ columns exist in protein group data
ibaq_cols = [col for col in ds.pg.data.columns
            if 'ibaq' in col.lower() and col.startswith('sample_')]

if ibaq_cols:
    # Plot first sample as example
    sample_col = ibaq_cols[0]
    ibaq_values = ds.pg.data[sample_col].dropna()

    if len(ibaq_values) > 0:
        # Create histogram with KDE
        fig, ax = plt.subplots(figsize=(10, 6))
        log_ibaq = np.log10(ibaq_values[ibaq_values > 0])
        log_ibaq.plot.hist(bins=50, alpha=0.6, ax=ax, label='Histogram')
        log_ibaq.plot.kde(ax=ax, label='KDE', linewidth=2)

        ax.set_xlabel('log10(iBAQ Intensity)')
        ax.set_ylabel('Frequency / Density')
        ax.set_title(f'iBAQ Distribution - {sample_col}')
        ax.legend()
        plt.tight_layout()

        # Save
        fig.savefig("./plots/ibaq_distribution.svg", format='svg', bbox_inches='tight')
```

**Output**:

- Histogram + kernel density estimate of iBAQ values
- Log-transformed for better visualization
- X-axis: log10(iBAQ intensity)
- Y-axis: Density or frequency

**Interpretation**:

- **Bimodal distribution**: Distinct protein abundance classes
- **Long tail**: High-abundance proteins (housekeeping, structural)
- **Narrow range**: Limited dynamic range, possible detection issues

## Missing Value Patterns

Visualize missing data patterns:

```python
import qpx
import matplotlib.pyplot as plt
import seaborn as sns

ds = qpx.Dataset("./output/")

# Get intensity columns
intensity_cols = [col for col in ds.feature.data.columns
                 if col.startswith('sample_')]

# Calculate missing value percentages
missing_pct = ds.feature.data[intensity_cols].isna().mean() * 100

# Create bar plot
fig, ax = plt.subplots(figsize=(12, 6))
missing_pct.plot.bar(ax=ax, color='coral')
ax.set_xlabel('Sample')
ax.set_ylabel('Missing Values (%)')
ax.set_title('Missing Value Percentage by Sample')
ax.axhline(y=30, color='r', linestyle='--', label='30% threshold')
plt.xticks(rotation=45, ha='right')
plt.legend()
plt.tight_layout()

# Save
fig.savefig("./plots/missing_values.svg", format='svg', bbox_inches='tight')
```

**Interpretation**:

- **Low missing values (<20%)**: Good data quality
- **Moderate (20-40%)**: Acceptable for most analyses
- **High (>40%)**: May require imputation or filtering

## Customization and Styling

### Publication-Quality Plots

Enhance plots for publication:

```python
import qpx
import matplotlib.pyplot as plt

# Set publication style
plt.style.use('seaborn-v0_8-paper')
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 16

ds = qpx.Dataset("./output/")

# Create plot with custom styling
fig = ds.identifications.plot()

# Save in multiple formats
fig.savefig("./plots/identifications.svg", format='svg', bbox_inches='tight', dpi=300)
fig.savefig("./plots/identifications.pdf", format='pdf', bbox_inches='tight', dpi=300)
fig.savefig("./plots/identifications.png", format='png', bbox_inches='tight', dpi=300)
```

### Multi-Panel Figures

Create comprehensive visualization panels:

```python
import qpx
import matplotlib.pyplot as plt

ds = qpx.Dataset("./output/")

# Create multi-panel figure
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Panel 1: Identifications
ax1 = axes[0, 0]
stats = {
    'Proteins': ds.psm.data['protein_accessions'].nunique(),
    'Peptides': ds.psm.data['sequence'].nunique(),
    'PSMs': ds.psm.count()
}
ax1.bar(stats.keys(), stats.values())
ax1.set_ylabel('Count')
ax1.set_title('A. Identifications')

# Panel 2: Runs
ax2 = axes[0, 1]
run_psms = ds.psm.data.groupby('run_file_name').size()
ax2.bar(range(len(run_psms)), run_psms.values)
ax2.set_xlabel('Run')
ax2.set_ylabel('PSM Count')
ax2.set_title('B. PSMs per Run')

# Panel 3: Modifications
ax3 = axes[1, 0]
if 'modifications' in ds.psm.data.columns:
    mod_counts = ds.psm.data['modifications'].value_counts().head(10)
    mod_counts.plot.barh(ax=ax3)
    ax3.set_xlabel('Count')
    ax3.set_title('C. Top 10 Modifications')

# Panel 4: Intensity distribution (if feature data available)
ax4 = axes[1, 1]
if hasattr(ds, 'feature') and ds.feature.count() > 0:
    intensity_cols = [col for col in ds.feature.data.columns
                     if col.startswith('sample_')][:5]
    for col in intensity_cols:
        intensities = ds.feature.data[col].dropna()
        if len(intensities) > 0:
            intensities[intensities > 0].apply('log10').plot.kde(ax=ax4, label=col)
    ax4.set_xlabel('log10(Intensity)')
    ax4.set_ylabel('Density')
    ax4.set_title('D. Intensity Distribution')
    ax4.legend()

plt.tight_layout()
fig.savefig("./plots/comprehensive_qc.svg", format='svg', bbox_inches='tight', dpi=300)
```

## General Plotting Tips

### Output Formats

- **SVG**: Recommended for publications (scalable, editable)
- **PDF**: Alternative vector format (portable)
- **PNG**: Raster format (use high DPI: 300+)

### Color Considerations

- Plots use color-blind friendly palettes by default
- Ensure sufficient contrast for grayscale printing
- Use ColorBrewer palettes for multi-category plots

### Size and Resolution

- Vector formats (SVG/PDF) scale without quality loss
- For raster formats, use DPI ≥ 300 for publications
- Standard figure sizes: single column (3.5"), double column (7")

### Best Practices

- **Label axes clearly**: Include units where applicable
- **Add titles**: Descriptive but concise
- **Use legends**: When plotting multiple series
- **Limit colors**: Use 5-10 distinct colors maximum
- **Save originals**: Keep vector formats for future editing

## Integration with Analysis Workflows

Combine statistics and visualization:

```python
import qpx
import matplotlib.pyplot as plt

def comprehensive_qc_workflow(dataset_path, output_dir):
    """Generate comprehensive QC report with plots."""
    ds = qpx.Dataset(dataset_path)

    # Generate all plots
    plots = {
        'identifications': ds.identifications.plot(),
        'runs': ds.runs.plot(),
        'modifications': ds.modifications.plot(),
        'qc': ds.qc.plot()
    }

    # Save plots
    for name, fig in plots.items():
        fig.savefig(f"{output_dir}/{name}.svg", format='svg', bbox_inches='tight')
        plt.close(fig)

    # Generate statistics summary
    with open(f"{output_dir}/statistics.txt", 'w') as f:
        f.write("QPX Quality Control Report\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Proteins: {ds.psm.data['protein_accessions'].nunique():,}\n")
        f.write(f"Peptides: {ds.psm.data['sequence'].nunique():,}\n")
        f.write(f"PSMs: {ds.psm.count():,}\n")
        f.write(f"Runs: {ds.psm.data['run_file_name'].nunique()}\n")

    print(f"QC report generated in: {output_dir}")

# Usage
comprehensive_qc_workflow("./output/", "./qc_report/")
```

## Advanced Customization

For advanced customization beyond the built-in views, you can access the underlying data directly:

```python
import qpx
import matplotlib.pyplot as plt
import pandas as pd

ds = qpx.Dataset("./output/")

# Access raw data
psm_df = ds.psm.data
feature_df = ds.feature.data
pg_df = ds.pg.data

# Create custom visualizations using matplotlib, seaborn, plotly, etc.
# Example: Custom scatter plot
fig, ax = plt.subplots(figsize=(10, 6))
ax.scatter(psm_df['rt'], psm_df['observed_mz'], alpha=0.5, s=10)
ax.set_xlabel('Retention Time (seconds)')
ax.set_ylabel('Observed m/z')
ax.set_title('PSM Distribution in RT-m/z Space')
plt.tight_layout()
fig.savefig("./plots/custom_scatter.svg", format='svg', bbox_inches='tight')
```

## Related Documentation

- [Statistics Guide](stats.md) - Generate numeric summaries
- [Dataset API Reference](../spec/dataset.md) - Core dataset structure
- [Views Specification](../spec/views.md) - Available dataset views
- [Quickstart Guide](../quickstart.md) - Get started with QPX
