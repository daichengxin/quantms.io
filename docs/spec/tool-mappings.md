# Tool Field Mappings

!!! info "Reference for converter developers"
    This page is a reference for converter developers. For schema definitions, see the individual view pages: [PSM](psm.md), [Feature](feature.md), [Protein Group](pg.md).

This page documents how fields from different proteomics search engines and tools map to the corresponding QPX fields. Each section covers a specific QPX view and lists the source column name in each supported tool.

A dash (`---`) indicates that the tool does not provide a direct mapping for that field.

## PSM Field Mappings

The PSM (Peptide Spectrum Match) view captures spectrum-level identification results. The table below shows how each QPX PSM field maps to columns in common proteomics tools.

| QPX Field | MaxQuant | DIA-NN | FragPipe | mzTab |
|---|---|---|---|---|
| `sequence` | Sequence | Stripped.Sequence | Peptide | sequence |
| `peptidoform` | Modified sequence | Modified.Sequence | Modified Peptide | opt_global_cv_MS:1000889_peptidoform_sequence |
| `charge` | Charge | Precursor.Charge | Charge | charge |
| `posterior_error_probability` | PEP | PEP | --- | opt_global_Posterior_Error_Probability_score |
| `is_decoy` | Reverse | --- | --- | opt_global_cv_MS:1002217_decoy_peptide |
| `calculated_mz` | --- | --- | Calculated M/Z | calc_mass_to_charge |
| `observed_mz` | m/z | --- | Observed M/Z | exp_mass_to_charge |
| `protein_accessions` | Proteins | Protein.Group | Protein | accession |
| `run_file_name` | Raw file | Run | Spectrum File | spectra_ref |
| `scan` | Scan number | --- | Scan | --- |
| `rt` | Retention time | RT | Retention | retention_time |

!!! note "PEP direction"
    `posterior_error_probability` is the probability that the PSM is incorrect — **lower values indicate higher confidence** (lower is better). All major tools (Percolator, MaxQuant) export PEP directly as P(incorrect). FragPipe exports PeptideProphet Probability (P(correct)), so PEP must be computed as `1 - probability`.

## Feature Field Mappings

The Feature view captures quantified peptide features with intensity data. Features aggregate information across scans and are the primary view for DIA and label-free quantification workflows.

| QPX Field | MaxQuant | DIA-NN | FragPipe | mzTab |
|---|---|---|---|---|
| `sequence` | Sequence | Stripped.Sequence | Peptide | sequence |
| `peptidoform` | Modified sequence | Modified.Sequence | Modified Peptide | opt_global_cv_MS:1000889_peptidoform_sequence |
| `charge` | Charge | Precursor.Charge | --- | charge |
| `is_decoy` | Reverse | --- | --- | opt_global_cv_MS:1002217_decoy_peptide |
| `calculated_mz` | --- | --- | Calculated M/Z | calc_mass_to_charge |
| `observed_mz` | m/z | --- | --- | exp_mass_to_charge |
| `rt` | Retention time | RT | --- | retention_time |
| `rt_start` | --- | RT.Start | --- | --- |
| `rt_stop` | --- | RT.Stop | --- | --- |
| `predicted_rt` | --- | Predicted.RT | --- | --- |
| `intensities` | Intensity | Precursor.Quantity | Intensity | Intensity |
| `pg_accessions` | Proteins | Protein.Group | --- | accession |
| `anchor_protein` | --- | --- | --- | --- |
| `pg_positions` | --- | --- | --- | --- |
| `run_file_name` | Raw file | Run | --- | --- |

## Protein Group Field Mappings

The Protein Group view captures protein-level quantification and inference results. The table below covers the three main supported tools; mzTab protein group mapping is handled through the PSM and Feature views.

| QPX Field | MaxQuant | DIA-NN | FragPipe |
|---|---|---|---|
| `pg_accessions` | Protein IDs | Protein.Group | Group + Indistinguishable Proteins |
| `pg_names` | Protein names | Protein.Names | --- |
| `gg_accessions` | Gene names | Genes | --- |
| `run_file_name` | combined | Run | --- |
| `peptide_counts` | Unique peptides | Unique.Stripped.Peptides | Unique Peptides |
| `feature_counts` | MS/MS count | Precursor.Quantity | Precursor Ions |
| `global_qvalue` | Q-value | Global.PG.Q.Value | --- |
| `intensities` | Intensity, LFQ intensity | PG.Quantity | --- |
| `additional_intensities` | LFQ intensity, iBAQ | PG.Normalised, PG.MaxLFQ | --- |
| `is_decoy` | Reverse | --- | --- |
| `contaminant` | Potential contaminant | --- | --- |

!!! note "Tool version considerations"
    Column names may vary between tool versions. The mappings above reflect the commonly used versions at the time of writing:

    - **MaxQuant**: v2.x (`evidence.txt`, `proteinGroups.txt`)
    - **DIA-NN**: v1.8+ (`report.tsv`, `pg_matrix.tsv`)
    - **FragPipe**: v20+ (combined output files)
    - **mzTab**: v1.0 specification columns

!!! tip "Adding support for a new tool"
    To add a converter for a new tool, implement a mapping from the tool's output columns to the QPX fields listed in the tables above. See the [PSM](psm.md), [Feature](feature.md), and [Protein Group](pg.md) schema pages for field type requirements and nullability constraints.
