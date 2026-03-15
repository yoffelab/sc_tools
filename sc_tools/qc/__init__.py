"""
sc_tools.qc: Quality control utilities for spatial and single-cell omics.

- metrics: calculate_qc_metrics, filter_cells, filter_genes, highly_variable_genes (scanpy)
- spatial: spatially_variable_genes (squidpy); % mt / % hb per spot via calculate_qc_metrics
- plots: qc_2x2_grid, qc_spatial_multipage, cross-sample comparison plots
- sample_qc: per-sample metrics, pass/fail classification, spot filtering
- report: HTML QC report generation
- doublet: Solo-based doublet scoring for single-cell resolution modalities
"""

from . import doublet, metrics, plots, report, sample_qc, spatial
from .doublet import score_doublets_solo
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
    qc_pct_mt_per_sample,
    qc_sample_comparison_bar,
    qc_sample_scatter_matrix,
    qc_sample_violin_grouped,
    qc_scatter_counts_genes,
    qc_spatial_multipage,
    qc_violin_metrics,
)
from .report import (
    generate_all_qc_reports,
    generate_post_celltyping_report,
    generate_post_filter_report,
    generate_post_integration_report,
    generate_pre_filter_report,
    generate_qc_report,
    generate_segmentation_qc_report,
)
from .sample_qc import (
    apply_qc_filter,
    classify_samples,
    compute_sample_metrics,
    filter_spots,
    save_pass_fail_lists,
)
from .spatial import spatially_variable_genes, spatially_variable_genes_per_library

__all__ = [
    "metrics",
    "spatial",
    "plots",
    "sample_qc",
    "report",
    "doublet",
    "score_doublets_solo",
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
    "qc_pct_mt_per_sample",
    "qc_sample_comparison_bar",
    "qc_sample_violin_grouped",
    "qc_sample_scatter_matrix",
    "filter_spots",
    "compute_sample_metrics",
    "classify_samples",
    "save_pass_fail_lists",
    "apply_qc_filter",
    "generate_qc_report",
    "generate_pre_filter_report",
    "generate_post_filter_report",
    "generate_post_integration_report",
    "generate_post_celltyping_report",
    "generate_segmentation_qc_report",
    "generate_all_qc_reports",
]
