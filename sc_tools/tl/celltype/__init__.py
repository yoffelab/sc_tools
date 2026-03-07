"""
sc_tools.tl.celltype: Automated cell-type annotation.

Provides a unified dispatcher (annotate_celltypes) and a curated-map
applier (apply_celltype_map) with pluggable backends.

Backends (Tier 1-2, no/optional soft deps)
-------------------------------------------
sctype        -- ScType marker scoring (zero soft deps)
custom_gates  -- Hierarchical threshold gating (zero soft deps)
celltypist    -- CellTypist probabilistic model (soft dep: celltypist)
ensemble      -- Voting over prior annotation columns (zero soft deps)

Backends (Tier 3, stubs pending full implementation)
-----------------------------------------------------
scarches      -- scArches reference mapping (stub)
geneformer    -- Geneformer foundation model (stub)
scgpt         -- scGPT foundation model (stub)
singler       -- SingleR R-based annotation (stub)

Usage
-----
>>> from sc_tools.tl.celltype import annotate_celltypes, apply_celltype_map
>>> annotate_celltypes(adata, method="sctype", marker_db=my_db)
>>> apply_celltype_map(adata, mapping="metadata/celltype_map.json")
"""

# Register backends by importing side-effect modules
from . import (  # noqa: F401
    _celltypist,
    _custom_gates,
    _ensemble,
    _geneformer,
    _scarches,
    _scgpt,
    _sctype,
    _singler,
)
from ._base import list_celltype_methods, register_celltype_backend
from .annotate import annotate_celltypes
from .apply import apply_celltype_map

__all__ = [
    "annotate_celltypes",
    "apply_celltype_map",
    "list_celltype_methods",
    "register_celltype_backend",
]
