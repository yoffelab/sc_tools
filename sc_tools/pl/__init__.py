"""
sc_tools.pl: Plotting utilities for spatial omics data.

Following scanpy's API pattern, this module provides plotting functions:
- spatial: Spatial plot wrappers
- heatmaps: Heatmap/clustermap utilities
- statistical: Statistical annotations (bars, asterisks)
- volcano: Volcano plot utilities
"""

# Import plotting functions
from . import spatial
from . import heatmaps
from . import statistical
from . import volcano
from . import save as save_figs
from .save import save_figure

__all__ = ['spatial', 'heatmaps', 'statistical', 'volcano', 'save_figs', 'save_figure']
