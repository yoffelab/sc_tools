"""
sc_tools.gr._nhood — nhood_enrichment and interaction_matrix per-ROI wrappers.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import scipy.stats as ss

from ._aggregate import apply_bh_correction, combine_pvalues, cross_roi_zscore, unify_matrices
from ._utils import iter_rois

try:
    import squidpy as sq

    SQUIDPY_AVAILABLE = True
except ImportError:
    SQUIDPY_AVAILABLE = False
    sq = None


def nhood_enrichment(
    adata: ad.AnnData,
    cluster_key: str,
    library_key: str = "library_id",
    n_perms: int = 1000,
    **kwargs,
) -> None:
    """
    Run sq.gr.nhood_enrichment per ROI and aggregate results across ROIs.

    Per-ROI results are stored in adata.uns['gr']['nhood_enrichment']['per_roi'].
    Aggregated statistics (zscore_mean, zscore_std, cross_roi_zscore,
    pval_stouffer, pval_fdr_bh) are stored in the 'aggregated' sub-dict.

    Parameters
    ----------
    adata
        Multi-ROI AnnData. Must have spatial_connectivities computed.
    cluster_key
        obs column with cell type / cluster categorical labels.
    library_key
        obs column identifying each ROI.
    n_perms
        Number of permutations for sq.gr.nhood_enrichment.
    **kwargs
        Additional keyword arguments for sq.gr.nhood_enrichment.
    """
    if not SQUIDPY_AVAILABLE:
        raise ImportError("squidpy is required for sc_tools.gr.nhood_enrichment")

    per_roi_zscores: list[np.ndarray] = []
    per_roi_counts: list[np.ndarray] = []
    per_roi_cats: list[list[str]] = []
    per_roi_results: dict = {}

    for roi_id, roi in iter_rois(adata, library_key=library_key, cluster_key=cluster_key):
        try:
            sq.gr.nhood_enrichment(roi, cluster_key=cluster_key, n_perms=n_perms, **kwargs)
            ne_key = f"{cluster_key}_nhood_enrichment"
            ne = roi.uns[ne_key]
            zscore = np.asarray(ne["zscore"], dtype=float)
            count = np.asarray(ne["count"], dtype=float)
            cats = list(roi.obs[cluster_key].cat.categories)

            per_roi_zscores.append(zscore)
            per_roi_counts.append(count)
            per_roi_cats.append(cats)
            per_roi_results[roi_id] = {
                "zscore": zscore,
                "count": count,
                "cats": cats,
            }
        except Exception as e:
            per_roi_results[roi_id] = {"error": str(e)}

    # Unify matrices
    valid_roi_ids = [k for k in per_roi_results if "error" not in per_roi_results[k]]
    valid_zscores = [per_roi_results[r]["zscore"] for r in valid_roi_ids]
    valid_counts = [per_roi_results[r]["count"] for r in valid_roi_ids]
    valid_cats = [per_roi_results[r]["cats"] for r in valid_roi_ids]

    if valid_zscores:
        unified_z, global_cats = unify_matrices(valid_zscores, valid_cats, fill_value=np.nan)
        unified_c, _ = unify_matrices(valid_counts, valid_cats, fill_value=0.0)

        zscore_mean = np.nanmean(unified_z, axis=0)
        zscore_std = np.nanstd(unified_z, axis=0)
        z_cross = cross_roi_zscore(unified_z)

        # Convert per-ROI z-scores to two-tailed p-values
        pval_arr = 2.0 * ss.norm.sf(np.abs(unified_z))
        pval_stouffer = combine_pvalues(pval_arr, method="stouffer")
        pval_fdr_bh = apply_bh_correction(pval_stouffer)

        aggregated = {
            "cats": global_cats,
            "zscore_mean": zscore_mean,
            "zscore_std": zscore_std,
            "cross_roi_zscore": z_cross,
            "pval_stouffer": pval_stouffer,
            "pval_fdr_bh": pval_fdr_bh,
            "count_sum": np.nansum(unified_c, axis=0),
        }
    else:
        aggregated = {}

    adata.uns.setdefault("gr", {})["nhood_enrichment"] = {
        "per_roi": per_roi_results,
        "aggregated": aggregated,
    }


def interaction_matrix(
    adata: ad.AnnData,
    cluster_key: str,
    library_key: str = "library_id",
    **kwargs,
) -> None:
    """
    Run sq.gr.interaction_matrix per ROI and aggregate.

    Absent cell types are filled with 0 (missing = zero edges).
    Results: count_sum (sum across ROIs), prop_mean (mean row-normalised).

    Parameters
    ----------
    adata
        Multi-ROI AnnData. Must have spatial_connectivities computed.
    cluster_key
        obs column with cell type / cluster labels.
    library_key
        obs column identifying each ROI.
    **kwargs
        Additional keyword arguments for sq.gr.interaction_matrix.
    """
    if not SQUIDPY_AVAILABLE:
        raise ImportError("squidpy is required for sc_tools.gr.interaction_matrix")

    per_roi_counts: list[np.ndarray] = []
    per_roi_cats: list[list[str]] = []

    for _roi_id, roi in iter_rois(adata, library_key=library_key, cluster_key=cluster_key):
        try:
            sq.gr.interaction_matrix(roi, cluster_key=cluster_key, normalized=False, **kwargs)
            im_key = f"{cluster_key}_interactions"
            count = np.asarray(roi.uns[im_key], dtype=float)
            cats = list(roi.obs[cluster_key].cat.categories)
            per_roi_counts.append(count)
            per_roi_cats.append(cats)
        except Exception:
            pass

    if not per_roi_counts:
        adata.uns.setdefault("gr", {})["interaction_matrix"] = {}
        return

    unified_c, global_cats = unify_matrices(per_roi_counts, per_roi_cats, fill_value=0.0)

    count_sum = np.nansum(unified_c, axis=0)

    # Row-normalise each ROI then average
    row_sums = unified_c.sum(axis=2, keepdims=True)
    row_sums[row_sums == 0] = np.nan
    prop_per_roi = unified_c / row_sums
    prop_mean = np.nanmean(prop_per_roi, axis=0)

    adata.uns.setdefault("gr", {})["interaction_matrix"] = {
        "cats": global_cats,
        "count_sum": count_sum,
        "prop_mean": prop_mean,
    }
