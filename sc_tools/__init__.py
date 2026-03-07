"""
sc-tools: Generic reusable utilities for spatial omics analysis.

This package provides modular, reusable tools following scanpy's API pattern:
- sc_tools.pp: Preprocessing (normalization, integration, clustering, modality recipes)
- sc_tools.pl: Plotting utilities (spatial, heatmaps, statistical annotations)
- sc_tools.tl: Analysis tools (statistical testing, colocalization, etc.)
- sc_tools.qc: QC metrics, filters, spatially variable genes, QC report plotting
- sc_tools.data: Data loading, preprocessing, and I/O
- sc_tools.memory: Memory management and GPU utilities

All code in this package should be generic and not project-specific.
"""

__version__ = "0.1.0"

# Import main modules following scanpy pattern
from . import bm, ingest, memory, pl, pp, qc, storage, tl, utils, validate

# Registry is optional (requires sqlalchemy)
try:
    from . import registry
except ImportError:
    pass

__all__ = [
    "pp",
    "pl",
    "tl",
    "qc",
    "utils",
    "memory",
    "bm",
    "ingest",
    "validate",
    "storage",
    "registry",
]
