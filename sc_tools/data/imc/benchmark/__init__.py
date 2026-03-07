"""IMC segmentation benchmark data infrastructure.

Submodules:
- catalog: Dataset discovery and ROI cataloging
- public: Public dataset download and standardization
- prepare: ROI validation, probability map generation, channel extraction
- config: BenchmarkConfig dataclass with YAML I/O
"""

from .catalog import (
    IMCDatasetEntry,
    ROIRecord,
    build_benchmark_catalog,
    discover_internal_datasets,
    discover_rois,
)
from .config import BenchmarkConfig
from .prepare import (
    extract_channel_by_name,
    extract_dna_channels,
    generate_probability_map,
    normalize_imc_intensity,
    validate_roi_files,
)
from .public import (
    PUBLIC_DATASETS,
    download_public_dataset,
    standardize_public_dataset,
)

__all__ = [
    "IMCDatasetEntry",
    "ROIRecord",
    "build_benchmark_catalog",
    "discover_internal_datasets",
    "discover_rois",
    "BenchmarkConfig",
    "extract_channel_by_name",
    "extract_dna_channels",
    "generate_probability_map",
    "normalize_imc_intensity",
    "validate_roi_files",
    "PUBLIC_DATASETS",
    "download_public_dataset",
    "standardize_public_dataset",
]
