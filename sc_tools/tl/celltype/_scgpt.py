"""
scGPT foundation model backend stub.

Install: ``pip install 'sc-tools[foundation]'``

This module is a placeholder; full implementation is deferred until a
scGPT-based workflow is required by an active project.
"""

from __future__ import annotations

from ._base import register_celltype_backend


class ScGPTBackend:
    """scGPT transformer-based cell typing (stub; not yet implemented)."""

    @staticmethod
    def run(adata, *, cluster_key="leiden", store_proba=False, **kwargs):
        raise NotImplementedError(
            "The 'scgpt' backend is not yet implemented. "
            "Use 'sctype', 'celltypist', or 'custom_gates' instead."
        )


register_celltype_backend("scgpt", ScGPTBackend)
