"""IO Gateway: tiered data loading and memory estimation.

All heavy imports (h5py, anndata, psutil) are deferred to method bodies
so that ``import sc_tools.io`` does not trigger them (CLI-08 lazy import rule).

Public API
----------
IOGateway      -- tiered read dispatch with memory guard
DataTier       -- T1_METADATA / T2_SUMMARY / T3_FULL enum
read_h5ad_metadata -- h5py-only metadata extraction
estimate_from_h5   -- pre-load memory estimation
"""

from __future__ import annotations

from sc_tools.io.estimate import estimate_from_h5
from sc_tools.io.gateway import DataTier, IOGateway
from sc_tools.io.metadata import read_h5ad_metadata

__all__ = [
    "IOGateway",
    "DataTier",
    "read_h5ad_metadata",
    "estimate_from_h5",
]
