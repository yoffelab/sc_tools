"""
scArches reference-mapping backend stub.

Install: ``pip install 'sc-tools[celltyping]'``

This module is a placeholder; full implementation is deferred until a
scArches-based workflow is required by an active project.
"""

from __future__ import annotations

from ._base import register_celltype_backend


class ScArchesBackend:
    """scArches reference atlas mapping (stub; not yet implemented)."""

    @staticmethod
    def run(adata, *, cluster_key="leiden", store_proba=False, **kwargs):
        raise NotImplementedError(
            "The 'scarches' backend is not yet implemented. "
            "Use 'sctype', 'celltypist', or 'custom_gates' instead."
        )


register_celltype_backend("scarches", ScArchesBackend)
