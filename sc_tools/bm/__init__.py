"""
sc_tools.bm: Benchmarking modules for pipeline decision-making.

Submodules:
- segmentation: Cell segmentation quality metrics and comparison
- integration: Batch correction quality metrics (scib-style)
- mask_io: Mask format adapters (Cellpose, StarDist, CellProfiler, TIFF)
- report: HTML report generation for benchmarking results
"""

from __future__ import annotations

from .integration import (
    compare_integrations,
    compute_composite_score,
    compute_integration_metrics,
    run_integration_benchmark,
)
from .mask_io import (
    load_cellpose_mask,
    load_cellprofiler_mask,
    load_deepcell_mask,
    load_mask,
    load_stardist_mask,
    load_tiff_mask,
)
from .segment import (
    run_cellpose,
    run_stardist,
)
from .segmentation import (
    compare_segmentations,
    compute_detection_metrics,
    compute_marker_quality,
    compute_morphology_metrics,
    compute_segmentation_accuracy,
    compute_size_distribution,
    compute_spatial_coherence,
    score_segmentation,
)

__all__ = [
    # Segmentation
    "compute_morphology_metrics",
    "compute_marker_quality",
    "compute_spatial_coherence",
    "compute_size_distribution",
    "compute_detection_metrics",
    "compute_segmentation_accuracy",
    "score_segmentation",
    "compare_segmentations",
    # Integration
    "compute_integration_metrics",
    "compute_composite_score",
    "compare_integrations",
    "run_integration_benchmark",
    # Segment runners
    "run_cellpose",
    "run_stardist",
    # Mask I/O
    "load_mask",
    "load_cellpose_mask",
    "load_stardist_mask",
    "load_cellprofiler_mask",
    "load_deepcell_mask",
    "load_tiff_mask",
]
