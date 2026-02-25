"""
sc_tools.tl: Analysis tools for spatial omics data.

Following scanpy's API pattern, this module provides analysis functions:
- testing: Statistical testing (Mann-Whitney, FDR correction)
- colocalization: Spatial colocalization analysis (correlation, Moran's I, enrichment)
- deconvolution: Deconvolution utilities
"""

# Import tool functions
from . import colocalization, deconvolution, testing
from . import io as io_save
from .score_signature import score_signature

__all__ = ["testing", "colocalization", "deconvolution", "io_save", "score_signature"]
