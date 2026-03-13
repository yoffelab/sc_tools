"""
sc_tools.gr._ligrec — ligrec per-ROI wrapper.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

from ._aggregate import apply_bh_correction, combine_pvalues
from ._utils import iter_rois

try:
    import squidpy as sq

    SQUIDPY_AVAILABLE = True
except ImportError:
    SQUIDPY_AVAILABLE = False
    sq = None


def ligrec(
    adata: ad.AnnData,
    cluster_key: str,
    library_key: str = "library_id",
    interactions=None,
    n_perms: int = 1000,
    **kwargs,
) -> None:
    """
    Run sq.gr.ligrec per ROI and combine results.

    LR pair means and p-values are aligned across ROIs (union of LR pairs and
    cluster pairs). P-values are combined with Fisher method and BH-corrected.

    Parameters
    ----------
    adata
        Multi-ROI AnnData.
    cluster_key
        obs column with cell type / cluster labels.
    library_key
        obs column identifying each ROI.
    interactions
        Pre-computed interactions DataFrame or resource name accepted by squidpy.
        If None and no internet connectivity, a RuntimeError is raised.
    n_perms
        Number of permutations.
    **kwargs
        Additional keyword arguments for sq.gr.ligrec.
    """
    if not SQUIDPY_AVAILABLE:
        raise ImportError("squidpy is required for sc_tools.gr.ligrec")

    if interactions is None:
        import warnings

        warnings.warn(
            "interactions=None: squidpy will attempt to download OmniPath. "
            "Pass interactions=<DataFrame> to avoid network dependency.",
            UserWarning,
            stacklevel=2,
        )

    per_roi_means: list[pd.DataFrame] = []
    per_roi_pvals: list[pd.DataFrame] = []

    for _roi_id, roi in iter_rois(adata, library_key=library_key, cluster_key=cluster_key):
        try:
            sq_kwargs: dict = {"cluster_key": cluster_key, "n_perms": n_perms}
            if interactions is not None:
                sq_kwargs["interactions"] = interactions
            sq_kwargs.update(kwargs)

            sq.gr.ligrec(roi, **sq_kwargs)
            lr_key = f"{cluster_key}_ligrec"
            lr = roi.uns[lr_key]
            means_df = lr["means"]
            pvals_df = lr["pvalues"]
            per_roi_means.append(means_df)
            per_roi_pvals.append(pvals_df)
        except Exception:
            pass

    if not per_roi_means:
        adata.uns.setdefault("gr", {})["ligrec"] = {}
        return

    # Align to union index (LR pairs) and union columns (cluster pairs)
    all_index = list(dict.fromkeys(idx for df in per_roi_means for idx in df.index.tolist()))
    all_cols = list(dict.fromkeys(col for df in per_roi_means for col in df.columns.tolist()))

    def _align(dfs: list[pd.DataFrame], fill: float) -> np.ndarray:
        """Stack aligned DataFrames into (n_roi, n_lr, n_cluster) array."""
        arrays = []
        for df in dfs:
            aligned = df.reindex(index=all_index, columns=all_cols, fill_value=fill)
            arrays.append(aligned.values)
        return np.stack(arrays, axis=0)

    means_arr = _align(per_roi_means, fill=0.0)
    pvals_arr = _align(per_roi_pvals, fill=1.0)

    means_combined = np.nanmean(means_arr, axis=0)

    # Combine p-values with Fisher method
    pvals_combined = combine_pvalues(pvals_arr, method="fisher")
    pvals_fdr = apply_bh_correction(pvals_combined)

    means_df_out = pd.DataFrame(means_combined, index=all_index, columns=all_cols)
    pvals_df_out = pd.DataFrame(pvals_combined, index=all_index, columns=all_cols)
    pvals_fdr_df = pd.DataFrame(pvals_fdr, index=all_index, columns=all_cols)

    adata.uns.setdefault("gr", {})["ligrec"] = {
        "means": means_df_out,
        "pvalues": pvals_df_out,
        "pvalues_fdr_bh": pvals_fdr_df,
    }
