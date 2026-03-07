"""
Geneformer foundation model backend stub.

Install: ``pip install 'sc-tools[foundation]'``

This module is a placeholder; full implementation is deferred until a
Geneformer-based workflow is required by an active project.
"""

from __future__ import annotations

from ._base import register_celltype_backend


class GeneformerBackend:
    """Geneformer transformer-based cell typing (stub; not yet implemented)."""

    @staticmethod
    def run(adata, *, cluster_key="leiden", store_proba=False, **kwargs):
        raise NotImplementedError(
            "The 'geneformer' backend is not yet implemented. "
            "Use 'sctype', 'celltypist', or 'custom_gates' instead."
        )


register_celltype_backend("geneformer", GeneformerBackend)
