"""QPX converter layer -- tool-specific adapters that produce QPX Parquet files.

Each converter:
    1. Creates its own DuckDB connection (NOT shared with Dataset)
    2. Loads tool output into DuckDB (TSV, mzTab, etc.)
    3. Transforms via SQL into QPX schema
    4. Pipes results into a Writer (FeatureWriter, PsmWriter, PgWriter, etc.)
"""

from qpx.converters.base import BaseConverter
from qpx.converters.diann.converter import DiaNNConverter
from qpx.converters.diann.feature_adapter import DiannFeatureAdapter
from qpx.converters.diann.pg_adapter import DiannPgAdapter
from qpx.converters.fragpipe.converter import FragPipeConverter
from qpx.converters.fragpipe.feature_adapter import FragPipeFeatureAdapter
from qpx.converters.fragpipe.pg_adapter import FragPipePgAdapter
from qpx.converters.fragpipe.psm_adapter import FragPipePsmAdapter
from qpx.converters.maxquant.converter import MaxQuantConverter
from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter
from qpx.converters.maxquant.pg_adapter import MaxQuantPgAdapter
from qpx.converters.maxquant.psm_adapter import MaxQuantPsmAdapter
from qpx.converters.mzidentml.psm_adapter import MzIdentMLPsmAdapter
from qpx.converters.openms.converter import OpenMSConverter
from qpx.converters.orchestrator import BaseOrchestrator, build_dataset_record
from qpx.converters.quantms.converter import QuantMSConverter
from qpx.converters.quantms.feature_adapter import QuantmsFeatureAdapter
from qpx.converters.quantms.pg_adapter import QuantmsPgAdapter
from qpx.converters.quantms.psm_adapter import QuantmsPsmAdapter
from qpx.converters.sdrf import SdrfConverter
from qpx.converters.spectronaut.converter import SpectronautConverter
from qpx.converters.spectronaut.feature_adapter import SpectronautFeatureAdapter
from qpx.converters.spectronaut.pg_adapter import SpectronautPgAdapter

__all__ = [
    # Base
    "BaseConverter",
    "BaseOrchestrator",
    "build_dataset_record",
    # SDRF
    "SdrfConverter",
    # QuantMS (mzTab)
    "QuantmsPsmAdapter",
    "QuantmsFeatureAdapter",
    "QuantmsPgAdapter",
    "QuantMSConverter",
    # DIA-NN
    "DiannFeatureAdapter",
    "DiannPgAdapter",
    "DiaNNConverter",
    # MaxQuant
    "MaxQuantPsmAdapter",
    "MaxQuantFeatureAdapter",
    "MaxQuantPgAdapter",
    "MaxQuantConverter",
    # FragPipe
    "FragPipePsmAdapter",
    "FragPipeFeatureAdapter",
    "FragPipePgAdapter",
    "FragPipeConverter",
    # mzIdentML
    "MzIdentMLPsmAdapter",
    # OpenMS
    "OpenMSConverter",
    # Spectronaut
    "SpectronautFeatureAdapter",
    "SpectronautPgAdapter",
    "SpectronautConverter",
]
