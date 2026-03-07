"""
sc_tools.bm: Benchmarking modules for pipeline decision-making.

Submodules:
- segmentation: Cell segmentation quality metrics and comparison
- integration: Batch correction quality metrics (scib-style)
- mask_io: Mask format adapters (Cellpose, StarDist, CellProfiler, TIFF)
- segment: Cellpose, StarDist, DeepCell runners
- report: HTML report generation for benchmarking results
- runner: Benchmark orchestration across ROIs
- analysis: Cross-dataset statistics and generalization
- slurm: SLURM job generation for HPC
- cli: Command-line interface
- deepcell_runner: DeepCell/Mesmer wrapper
- strategy_dna: Strategy 2 (DNA-only pipeline)
- strategy_hf: Strategy 3 (HuggingFace pretrained models)
- strategy_vit: Strategy 4 (SegFormer trainable)
- postprocess: Shared post-processing utilities
"""

from __future__ import annotations

from .integration import (
    compare_integrations,
    compute_composite_score,
    compute_integration_metrics,
    run_full_integration_workflow,
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
    run_all_strategy1,
    run_cellpose,
    run_deepcell,
    run_stardist,
)
from .segmentation import (
    compare_segmentations,
    compute_boundary_metrics,
    compute_cell_type_preservation,
    compute_detection_metrics,
    compute_marker_quality,
    compute_morphology_metrics,
    compute_panoptic_quality,
    compute_segmentation_accuracy,
    compute_size_distribution,
    compute_spatial_coherence,
    score_segmentation,
)

__all__ = [
    # Segmentation metrics
    "compute_morphology_metrics",
    "compute_marker_quality",
    "compute_spatial_coherence",
    "compute_size_distribution",
    "compute_detection_metrics",
    "compute_segmentation_accuracy",
    "compute_panoptic_quality",
    "compute_boundary_metrics",
    "compute_cell_type_preservation",
    "score_segmentation",
    "compare_segmentations",
    # Integration
    "compute_integration_metrics",
    "compute_composite_score",
    "compare_integrations",
    "run_integration_benchmark",
    "run_full_integration_workflow",
    # Segment runners
    "run_cellpose",
    "run_stardist",
    "run_deepcell",
    "run_all_strategy1",
    # Mask I/O
    "load_mask",
    "load_cellpose_mask",
    "load_stardist_mask",
    "load_cellprofiler_mask",
    "load_deepcell_mask",
    "load_tiff_mask",
]
