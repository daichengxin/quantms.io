# Peptidoform

A **peptidoform** is a specific molecular form of a peptide: the amino acid sequence together with all of its modifications. In QPX, peptidoforms are encoded using the [ProForma notation](https://github.com/HUPO-PSI/ProForma), which provides a standardized, human-readable string representation for modified peptide sequences.

The `peptidoform` field appears alongside the bare `sequence` field (amino acids only) in every identification- and quantification-level view. While `sequence` gives the unmodified backbone, `peptidoform` captures the complete molecular identity needed to distinguish different modified forms of the same peptide.

## ProForma notation

ProForma is a HUPO-PSI standard for encoding peptide sequences with modifications in a single string. Modifications are placed in square brackets immediately after the residue they modify (or before the sequence for N-terminal modifications).

### Examples

**Basic modification by name**

```text
PEPTIDM[Oxidation]K
```

Methionine at position 7 carries an Oxidation modification.

**Modification by UNIMOD accession (recommended)**

```text
PEPTIDM[UNIMOD:35]K
```

Same as above, but using the UNIMOD accession number. This is the recommended approach because accessions are unambiguous.

**N-terminal modification**

```text
[Acetyl]-PEPTIDMK
```

Acetylation on the N-terminus of the peptide.

**Multiple modifications**

```text
[Acetyl]-PEPT[Phospho]IDM[Oxidation]K
```

N-terminal acetylation, phosphorylation on threonine at position 4, and oxidation on methionine at position 7.

**TMT-labeled peptide**

```text
[TMT6plex]-PEPTIDM[Oxidation]K[TMT6plex]
```

TMT6plex isobaric tag on the N-terminus and on the C-terminal lysine, with an additional oxidation on methionine.

!!! tip
    When a modification does not have a UNIMOD entry, use the CHEMMOD notation to specify the mass shift in Daltons: `PEPTIDM[CHEMMOD:+15.9949]K`.

## Where peptidoform is used

The `peptidoform` field is present in the following QPX views:

| View                  | Field name     | Description                                                    |
| --------------------- | -------------- | -------------------------------------------------------------- |
| PSM (`psm_file`)      | `peptidoform`  | The identified peptidoform for this spectrum match              |
| Feature (`feature_file`) | `peptidoform` | The quantified peptidoform for this feature                   |
| Peptide (`peptide_file`) | `peptidoform` | The peptidoform summarized across samples                     |

In all three views, the `peptidoform` string serves as a key component of the record identity alongside `charge` and `run_file_name`.

## Relationship to the modifications field

The `peptidoform` string encodes modifications inline, which is compact and human-readable. However, some use cases require structured access to modification details -- for example, querying all phosphorylation sites with their localization probabilities.

For this purpose, QPX provides a separate `modifications` field that represents the same information as a structured array of records. See the [Modifications](modifications.md) page for the full struct definition and examples.

!!! warning
    The `peptidoform` string and the `modifications` struct should be consistent with each other. If both are present, consumers should treat the `peptidoform` as the canonical representation and the `modifications` field as supplementary detail.

## Further reading

- [ProForma specification (HUPO-PSI)](https://github.com/HUPO-PSI/ProForma)
- [Modifications](modifications.md) -- structured modification representation with localization scores
- [QPX Format Overview](index.md) -- full list of views and concepts
