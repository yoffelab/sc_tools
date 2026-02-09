"""
sc_tools.tl: Analysis tools for spatial omics data.

Following scanpy's API pattern, this module provides analysis functions:
- testing: Statistical testing (Mann-Whitney, FDR correction)
- colocalization: Spatial colocalization analysis (correlation, Moran's I, enrichment)
- deconvolution: Deconvolution utilities
"""

# Import tool functions
from . import testing
from . import colocalization
from . import deconvolution
from . import io as io_save

__all__ = ['testing', 'colocalization', 'deconvolution', 'io_save']
