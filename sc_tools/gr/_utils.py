"""
Utility helpers for sc_tools.gr — per-ROI iteration and category alignment.
"""

from __future__ import annotations

import warnings
from collections.abc import Generator

import anndata as ad


def iter_rois(
    adata: ad.AnnData,
    library_key: str = "library_id",
    cluster_key: str | None = None,
    min_cells_per_type: int = 5,
) -> Generator[tuple[str, ad.AnnData], None, None]:
    """
    Yield (roi_id, roi_adata) for each ROI.

    Removes unused categories from cluster_key so that squidpy functions
    only see cell types actually present in the ROI.

    After removing unused categories, checks cell counts per category.
    If any category has fewer than min_cells_per_type cells, emits a
    UserWarning naming the ROI and the underpopulated category. The ROI
    is still yielded — callers decide whether to skip.

    Parameters
    ----------
    adata
        Full multi-ROI AnnData.
    library_key
        obs column containing ROI identifiers.
    cluster_key
        obs column with cell type / cluster labels (categorical).
        When provided, unused categories are removed per ROI.
    min_cells_per_type
        Minimum cells per category. Categories below this threshold trigger
        a UserWarning. Default: 5.
    """
    for roi_id in adata.obs[library_key].unique():
        roi = adata[adata.obs[library_key] == roi_id].copy()
        if cluster_key is not None and cluster_key in roi.obs.columns:
            roi.obs[cluster_key] = (
                roi.obs[cluster_key].astype("category").cat.remove_unused_categories()
            )
            # Warn for underpopulated categories
            counts = roi.obs[cluster_key].value_counts()
            for cat, count in counts.items():
                if count < min_cells_per_type:
                    warnings.warn(
                        f"ROI '{roi_id}': category '{cat}' has {count} cells, "
                        f"which is below min_cells_per_type={min_cells_per_type}.",
                        UserWarning,
                        stacklevel=2,
                    )
        yield roi_id, roi
