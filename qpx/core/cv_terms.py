"""PSI-MS controlled vocabulary term constants.

Centralizes CV accession strings used across QPX converters and core modules,
avoiding duplication and providing a single source of truth.
"""

# ---------------------------------------------------------------------------
# Cross-linking (mzIdentML 1.3 XL-MS)
# ---------------------------------------------------------------------------
CV_XL_DONOR = "MS:1002509"  # crosslink donor
CV_XL_ACCEPTOR = "MS:1002510"  # crosslink acceptor / receiver
CV_CROSSLINK_SII = "MS:1002511"  # crosslink spectrum identification item
CV_LOOPLINK_SII = "MS:1003329"  # looplink spectrum identification item
CV_NONCOVALENT_SII = "MS:1003331"  # noncovalently associated peptides SII

# ---------------------------------------------------------------------------
# Identification metadata
# ---------------------------------------------------------------------------
CV_PEPTIDE_UNIQUE = "MS:1001363"  # peptide unique to one protein
CV_PEPTIDE_SHARED = "MS:1001175"  # peptide shared in multiple proteins
CV_PEPTIDE_PASSES_THRESHOLD = "MS:1002500"  # peptide passes threshold
CV_PEPTIDOFORM_SEQUENCE = "MS:1000889"  # peptidoform sequence
CV_DECOY_PEPTIDE = "MS:1002217"  # decoy peptide

# ---------------------------------------------------------------------------
# Score direction indicators
# ---------------------------------------------------------------------------
CV_HIGHER_BETTER = "MS:1002108"  # higher score better
CV_LOWER_BETTER = "MS:1002109"  # lower score better

# ---------------------------------------------------------------------------
# Score parent terms (OBO hierarchy)
# ---------------------------------------------------------------------------
SCORE_PARENT_IDS = frozenset(
    {
        "MS:1001143",  # PSM-level search engine specific statistic
        "MS:1001153",  # search engine specific score
        "MS:1002347",  # PSM-level identification statistic
        "MS:1002363",  # search engine specific score for proteins
        "MS:1002368",  # search engine specific score for protein groups
        "MS:1002701",  # PSM-level result list statistic
        "MS:1002906",  # search engine specific score for proteoforms
    }
)

# ---------------------------------------------------------------------------
# Scores -- specific score accessions
# ---------------------------------------------------------------------------
CV_MASCOT_EXPECTATION = "MS:1001172"  # Mascot expectation value
CV_PSM_QVALUE = "MS:1002354"  # PSM-level q-value
CV_PSM_GLOBAL_FDR = "MS:1002350"  # PSM-level global FDR
CV_PROTEIN_GLOBAL_FDR = "MS:1001214"  # protein-level global FDR
CV_XL_PSM_GLOBAL_FDR = "MS:1003337"  # crosslinked PSM-level global FDR
CV_XL_PAIR_GLOBAL_FDR = "MS:1003338"  # peptide-pair sequence-level global FDR
CV_PEP = "MS:1001493"  # posterior error probability (PEP)
CV_PEP_GLOBAL = "MS:1002352"  # PSM-level global FDR (global PEP)

# ---------------------------------------------------------------------------
# Site localization scores
# ---------------------------------------------------------------------------
CV_PHOSPHORS_SITE_PROB = "MS:1001969"  # phosphoRS: site probability
CV_PTMRS_SITE_PROB = "MS:1002536"  # PTM-score (D-score / Ascore)
CV_ASCORE = "MS:1001985"  # Ascore

# Set of accessions that represent per-site localization scores.
# When found on a SpectrumIdentificationItem these should be mapped
# to the modification positions in the PSM record.
SITE_LOCALIZATION_ACCESSIONS = frozenset(
    {
        CV_PHOSPHORS_SITE_PROB,
        CV_PTMRS_SITE_PROB,
        CV_ASCORE,
    }
)

# ---------------------------------------------------------------------------
# Spectrum metadata
# ---------------------------------------------------------------------------
CV_SCAN_START_TIME = "MS:1000016"  # scan start time
CV_MISSED_CLEAVAGES = "MS:1003044"  # number of missed cleavages

# ---------------------------------------------------------------------------
# Derived sets (convenience groupings for converter logic)
# ---------------------------------------------------------------------------
SKIP_SCORE_ACCESSIONS = frozenset(
    {
        CV_PEPTIDE_UNIQUE,
        CV_PEPTIDE_SHARED,
        CV_PEPTIDE_PASSES_THRESHOLD,
        CV_CROSSLINK_SII,
        CV_LOOPLINK_SII,
        CV_NONCOVALENT_SII,
    }
)

LOWER_IS_BETTER_ACCESSIONS = frozenset(
    {
        CV_MASCOT_EXPECTATION,
        CV_PSM_QVALUE,
        CV_XL_PSM_GLOBAL_FDR,
        CV_XL_PAIR_GLOBAL_FDR,
    }
)
