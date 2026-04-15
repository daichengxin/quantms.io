"""Central constants for QPX structure and view names.

Use these instead of magic strings like "psm", "feature", "pg" to reduce
drift and typos. Align with Dataset._STRUCTURE_REGISTRY keys.
"""

# Structure / view names for ontology, provenance, and output
PSM = "psm"
FEATURE = "feature"
PG = "pg"
SAMPLE = "sample"
RUN = "run"
DATASET = "dataset"
ONTOLOGY = "ontology"
PROVENANCE = "provenance"
PEPMAP = "pepmap"
MZ = "mz"

# Common structure sets for converters
QUANT_STRUCTURES = (PSM, FEATURE, PG)
SAMPLE_RUN_STRUCTURES = (SAMPLE, RUN)
METADATA_STRUCTURES = (SAMPLE, RUN, ONTOLOGY)
