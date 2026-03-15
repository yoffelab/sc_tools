"""
sc_tools.gr._centrality — centrality_scores per-ROI wrapper.
"""

from __future__ import annotations

import warnings

import anndata as ad
import numpy as np
import pandas as pd

from ._aggregate import unify_dataframes
from ._utils import iter_rois

try:
    import squidpy as sq

    SQUIDPY_AVAILABLE = True
except ImportError:
    SQUIDPY_AVAILABLE = False
    sq = None


def centrality_scores(
    adata: ad.AnnData,
    cluster_key: str,
    library_key: str = "library_id",
    **kwargs,
) -> None:
    """
    Run sq.gr.centrality_scores per ROI and aggregate.

    Always calls cat.remove_unused_categories() before each ROI run.
    Results stored in adata.uns['gr']['centrality_scores'] with 'mean' and 'std'.

    Parameters
    ----------
    adata
        Multi-ROI AnnData. Must have spatial_connectivities computed.
    cluster_key
        obs column with cell type / cluster labels.
    library_key
        obs column identifying each ROI.
    **kwargs
        Additional keyword arguments for sq.gr.centrality_scores.
    """
    if not SQUIDPY_AVAILABLE:
        raise ImportError("squidpy is required for sc_tools.gr.centrality_scores")

    per_roi_dfs: list[pd.DataFrame] = []

    for roi_id, roi in iter_rois(adata, library_key=library_key, cluster_key=cluster_key):
        try:
            sq.gr.centrality_scores(roi, cluster_key=cluster_key, **kwargs)
            cs_key = f"{cluster_key}_centrality_scores"
            cs_df = roi.uns[cs_key]
            if not isinstance(cs_df, pd.DataFrame):
                continue
            per_roi_dfs.append(cs_df)
        except Exception as exc:
            warnings.warn(
                f"ROI '{roi_id}': centrality_scores failed: {exc}",
                UserWarning,
                stacklevel=2,
            )

    if not per_roi_dfs:
        adata.uns.setdefault("gr", {})["centrality_scores"] = {}
        return

    aligned = unify_dataframes(per_roi_dfs, fill_value=np.nan)
    stacked = np.stack([df.values for df in aligned], axis=0)  # (n_roi, n_cats, n_metrics)
    mean_vals = np.nanmean(stacked, axis=0)
    std_vals = np.nanstd(stacked, axis=0)

    cols = aligned[0].columns
    idx = aligned[0].index

    mean_df = pd.DataFrame(mean_vals, index=idx, columns=cols)
    std_df = pd.DataFrame(std_vals, index=idx, columns=cols)

    adata.uns.setdefault("gr", {})["centrality_scores"] = {
        "mean": mean_df,
        "std": std_df,
    }
