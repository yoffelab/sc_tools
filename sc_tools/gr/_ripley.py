"""
sc_tools.gr._ripley — Ripley statistics per-ROI wrapper.

CRITICAL: ripley must NEVER run on full concatenated adata.
"""

from __future__ import annotations

import warnings

import anndata as ad
import numpy as np

from ._aggregate import apply_bh_correction, combine_pvalues
from ._utils import iter_rois

try:
    import squidpy as sq

    SQUIDPY_AVAILABLE = True
except ImportError:
    SQUIDPY_AVAILABLE = False
    sq = None


def ripley(
    adata: ad.AnnData,
    cluster_key: str,
    library_key: str = "library_id",
    mode: str = "L",
    **kwargs,
) -> None:
    """
    Run sq.gr.ripley per ROI only (never on full concatenated adata).

    Always calls cat.remove_unused_categories() before each run.
    Per-ROI results are stored in adata.uns['gr']['ripley']['per_roi'].
    Where squidpy returns per-cluster p-values, these are combined across ROIs
    with the Stouffer method and BH-corrected, stored in 'combined'. Ripley
    statistics themselves are not averaged across ROIs because they are
    spatially bounded and distance-dependent; the p-value combination is the
    appropriate cross-ROI summary.

    Parameters
    ----------
    adata
        Multi-ROI AnnData.
    cluster_key
        obs column with cell type / cluster labels.
    library_key
        obs column identifying each ROI.
    mode
        Ripley statistic mode ('F', 'G', 'L').
    **kwargs
        Additional keyword arguments for sq.gr.ripley.
    """
    if not SQUIDPY_AVAILABLE:
        raise ImportError("squidpy is required for sc_tools.gr.ripley")

    per_roi: dict = {}
    per_cluster_pvals: dict[str, list[np.ndarray]] = {}

    for roi_id, roi in iter_rois(adata, library_key=library_key, cluster_key=cluster_key):
        try:
            sq.gr.ripley(roi, cluster_key=cluster_key, mode=mode, **kwargs)
            rip_key = f"{cluster_key}_ripley_{mode}"
            rip = roi.uns[rip_key]
            per_roi[roi_id] = rip

            # pvalues shape: (n_clusters, n_bins)
            pvalues = np.asarray(rip.get("pvalues", []), dtype=float)
            obs_df = rip.get(f"{mode}_stat")
            if obs_df is not None and pvalues.ndim == 2:
                clusters = list(obs_df.columns) if hasattr(obs_df, "columns") else []
                for c_idx, cluster in enumerate(clusters):
                    pvals = pvalues[c_idx]
                    per_cluster_pvals.setdefault(cluster, []).append(pvals)
        except Exception as exc:
            warnings.warn(
                f"ROI '{roi_id}': ripley failed: {exc}",
                UserWarning,
                stacklevel=2,
            )
            per_roi[roi_id] = {"error": str(exc)}

    # Combine p-values per cluster x distance bin
    combined: dict = {}
    for cluster, pval_list in per_cluster_pvals.items():
        if not pval_list:
            continue
        # Pad to same length
        max_len = max(len(p) for p in pval_list)
        padded = np.full((len(pval_list), max_len), 1.0)
        for i, p in enumerate(pval_list):
            padded[i, : len(p)] = p

        pval_comb = combine_pvalues(padded.reshape(len(pval_list), 1, max_len), method="stouffer")[
            0
        ]
        pval_fdr = apply_bh_correction(pval_comb)
        combined[cluster] = {
            "pval_stouffer": pval_comb,
            "pval_fdr_bh": pval_fdr,
        }

    adata.uns.setdefault("gr", {})["ripley"] = {
        "per_roi": per_roi,
        "combined": combined,
    }
