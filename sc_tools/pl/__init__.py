"""
sc_tools.pl: Plotting utilities for spatial omics data.

Following scanpy's API pattern, this module provides plotting functions:
- spatial: Spatial plot wrappers
- heatmaps: Heatmap/clustermap utilities
- statistical: Statistical annotations (bars, asterisks)
- volcano: Volcano plot utilities
- QC plots: qc_2x2_grid, qc_spatial_multipage (from sc_tools.qc.plots)
"""

# Import plotting functions
from . import heatmaps, spatial, statistical, volcano
from . import save as save_figs
from .save import save_figure

# QC plots (Phase 1): re-export from sc_tools.qc.plots for st.pl.qc_* usage
from sc_tools.qc.plots import (
    plot_highly_variable_genes,
    plot_spatially_variable_genes,
    qc_2x2_grid,
    qc_2x4_pre_post,
    qc_scatter_counts_genes,
    qc_spatial_multipage,
    qc_violin_metrics,
)

__all__ = [
    "spatial",
    "heatmaps",
    "statistical",
    "volcano",
    "save_figs",
    "save_figure",
    "qc_2x2_grid",
    "qc_2x4_pre_post",
    "qc_spatial_multipage",
    "qc_violin_metrics",
    "qc_scatter_counts_genes",
    "plot_highly_variable_genes",
    "plot_spatially_variable_genes",
]
