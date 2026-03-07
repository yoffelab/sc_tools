"""
SingleR reference-based annotation backend stub.

Install: Requires R + SingleR via rpy2.

This module is a placeholder; full implementation is deferred until a
SingleR-based workflow is required by an active project.
"""

from __future__ import annotations

from ._base import register_celltype_backend


class SingleRBackend:
    """SingleR R-based reference annotation (stub; not yet implemented)."""

    @staticmethod
    def run(adata, *, cluster_key="leiden", store_proba=False, **kwargs):
        raise NotImplementedError(
            "The 'singler' backend is not yet implemented. "
            "Use 'sctype', 'celltypist', or 'custom_gates' instead."
        )


register_celltype_backend("singler", SingleRBackend)
