"""
sc_tools.tl: Analysis tools for spatial omics data.

Following scanpy's API pattern, this module provides analysis functions:
- testing: Statistical testing (Mann-Whitney, FDR correction)
- colocalization: Spatial colocalization analysis (correlation, Moran's I, enrichment)
- deconvolution: Cell-type deconvolution (cell2location, tangram, destvi)
"""

# Import tool functions
from . import colocalization, testing
from . import io as io_save
from .deconvolution import (
    deconvolution,
    extract_reference_profiles,
    select_signature_genes,
)
from .score_signature import score_signature

__all__ = [
    "testing",
    "colocalization",
    "io_save",
    "score_signature",
    "deconvolution",
    "extract_reference_profiles",
    "select_signature_genes",
]
