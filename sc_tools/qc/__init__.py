"""
sc_tools.qc: Quality control utilities for spatial and single-cell omics.

- metrics: calculate_qc_metrics, filter_cells, filter_genes, highly_variable_genes (scanpy)
- spatial: spatially_variable_genes (squidpy); % mt / % hb per spot via calculate_qc_metrics
- plots: qc_2x2_grid, qc_spatial_multipage (2x2 metric grid; multipage spatial total_count, log1p, %mt)
"""

from . import metrics, plots, spatial
from .metrics import (
    calculate_qc_metrics,
    filter_cells,
    filter_genes,
    highly_variable_genes,
)
from .plots import (
    plot_highly_variable_genes,
    plot_spatially_variable_genes,
    qc_2x2_grid,
    qc_2x4_pre_post,
    qc_scatter_counts_genes,
    qc_spatial_multipage,
    qc_violin_metrics,
)
from .spatial import spatially_variable_genes, spatially_variable_genes_per_library

__all__ = [
    "metrics",
    "spatial",
    "plots",
    "calculate_qc_metrics",
    "filter_cells",
    "filter_genes",
    "highly_variable_genes",
    "spatially_variable_genes",
    "spatially_variable_genes_per_library",
    "qc_2x2_grid",
    "qc_2x4_pre_post",
    "qc_spatial_multipage",
    "qc_violin_metrics",
    "qc_scatter_counts_genes",
    "plot_highly_variable_genes",
    "plot_spatially_variable_genes",
]
