
# ============================================================================
# chembl_curator/__init__.py
# ============================================================================

"""ChEMBL Curator: A package for curating ChEMBL bioactivity data."""

__version__ = "0.1.0"
__author__ = "Your Name"

from .curator import ChEMBLCurator
from .downloader import ChEMBLDownloader
from .filters import ActivityFilter, CompoundFilter
from .config import CurationConfig

__all__ = [
    "ChEMBLCurator",
    "ChEMBLDownloader", 
    "ActivityFilter",
    "CompoundFilter",
    "CurationConfig",
]

