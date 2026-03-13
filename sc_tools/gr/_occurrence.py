"""
sc_tools.gr._occurrence — co_occurrence per-ROI wrapper.

CRITICAL: co_occurrence must NEVER run on full concatenated adata.
Coordinates from different ROIs must not mix.
"""

from __future__ import annotations

import anndata as ad
import numpy as np

from ._aggregate import cross_roi_zscore
from ._utils import iter_rois

try:
    import squidpy as sq

    SQUIDPY_AVAILABLE = True
except ImportError:
    SQUIDPY_AVAILABLE = False
    sq = None


def co_occurrence(
    adata: ad.AnnData,
    cluster_key: str,
    library_key: str = "library_id",
    **kwargs,
) -> None:
    """
    Run sq.gr.co_occurrence per ROI only (never on full concatenated adata).

    Results stored in adata.uns['gr']['co_occurrence']['per_roi'] and
    aggregated occ_mean, occ_std, occ_z in 'aggregated'.

    Parameters
    ----------
    adata
        Multi-ROI AnnData.
    cluster_key
        obs column with cell type / cluster labels.
    library_key
        obs column identifying each ROI.
    **kwargs
        Additional keyword arguments for sq.gr.co_occurrence.
    """
    if not SQUIDPY_AVAILABLE:
        raise ImportError("squidpy is required for sc_tools.gr.co_occurrence")

    per_roi: dict = {}
    occ_arrays: list[np.ndarray] = []
    cat_lists: list[list[str]] = []
    interval_ref: np.ndarray | None = None

    for roi_id, roi in iter_rois(adata, library_key=library_key, cluster_key=cluster_key):
        try:
            sq.gr.co_occurrence(roi, cluster_key=cluster_key, **kwargs)
            co_key = f"{cluster_key}_co_occurrence"
            co = roi.uns[co_key]
            occ = np.asarray(co["occ"], dtype=float)  # (N, N, D)
            interval = np.asarray(co["interval"], dtype=float)
            cats = list(roi.obs[cluster_key].cat.categories)

            per_roi[roi_id] = {
                "occ": occ,
                "interval": interval,
                "cats": cats,
            }
            occ_arrays.append(occ)
            cat_lists.append(cats)
            if interval_ref is None:
                interval_ref = interval
        except Exception as e:
            per_roi[roi_id] = {"error": str(e)}

    # Aggregate if we have results
    if occ_arrays:
        # Build global category list
        seen: dict[str, None] = {}
        for cats in cat_lists:
            for c in cats:
                seen[c] = None
        global_cats = list(seen.keys())
        cat_to_idx = {c: i for i, c in enumerate(global_cats)}
        n_global = len(global_cats)
        n_dist = occ_arrays[0].shape[2] if occ_arrays[0].ndim == 3 else 1
        n_roi = len(occ_arrays)

        unified = np.full((n_roi, n_global, n_global, n_dist), np.nan, dtype=float)
        for r, (occ, cats) in enumerate(zip(occ_arrays, cat_lists, strict=False)):
            for i, ci in enumerate(cats):
                gi = cat_to_idx[ci]
                for j, cj in enumerate(cats):
                    gj = cat_to_idx[cj]
                    unified[r, gi, gj, :] = occ[i, j, :]

        occ_mean = np.nanmean(unified, axis=0)
        occ_std = np.nanstd(unified, axis=0)

        # z-score per distance bin: treat (n_roi, N, N) for each bin
        occ_z = np.full_like(unified, np.nan)
        for d in range(n_dist):
            occ_z[:, :, :, d] = cross_roi_zscore(unified[:, :, :, d])

        aggregated = {
            "cats": global_cats,
            "interval": interval_ref,
            "occ_mean": occ_mean,
            "occ_std": occ_std,
            "occ_z": occ_z,
        }
    else:
        aggregated = {}

    adata.uns.setdefault("gr", {})["co_occurrence"] = {
        "per_roi": per_roi,
        "aggregated": aggregated,
    }
