"""
sc_tools.gr._autocorr — spatial_autocorr per-ROI wrapper.
"""

from __future__ import annotations

import warnings

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


def spatial_autocorr(
    adata: ad.AnnData,
    library_key: str = "library_id",
    mode: str = "moran",
    genes: list[str] | None = None,
    n_perms: int = 1000,
    n_jobs: int = 1,
    **kwargs,
) -> None:
    """
    Run sq.gr.spatial_autocorr per ROI and aggregate.

    Moran I is combined with Fisher r-to-z transform.
    P-values combined with Stouffer weighted by sqrt(n_cells_in_roi) and BH-corrected.

    Only mode='moran' is supported for cross-ROI aggregation. Geary C is bounded
    [0, 2] and cannot be aggregated with Fisher r-to-z (which requires [-1, 1]).

    Parameters
    ----------
    adata
        Multi-ROI AnnData. Must have spatial_connectivities computed.
    library_key
        obs column identifying each ROI.
    mode
        Must be "moran". "geary" raises NotImplementedError.
    genes
        Genes to test. Defaults to all var_names.
    n_perms
        Number of permutations.
    n_jobs
        Parallel jobs.
    **kwargs
        Additional keyword arguments for sq.gr.spatial_autocorr.
    """
    if not SQUIDPY_AVAILABLE:
        raise ImportError("squidpy is required for sc_tools.gr.spatial_autocorr")

    if mode == "geary":
        raise NotImplementedError(
            "Fisher r-to-z aggregation is not valid for Geary C (bounded [0, 2], not [-1, 1]). "
            "Use mode='moran' for cross-ROI aggregation."
        )

    per_roi_results: list[pd.DataFrame] = []
    per_roi_n_cells: list[int] = []

    for roi_id, roi in iter_rois(adata, library_key=library_key):
        try:
            sq.gr.spatial_autocorr(
                roi,
                mode=mode,
                genes=genes,
                n_perms=n_perms,
                n_jobs=n_jobs,
                **kwargs,
            )
            key = "moranI" if mode == "moran" else "gearyC"
            result_df = roi.uns[key]
            per_roi_results.append(result_df)
            per_roi_n_cells.append(roi.n_obs)
        except Exception as exc:
            warnings.warn(
                f"ROI '{roi_id}': spatial_autocorr failed: {exc}",
                UserWarning,
                stacklevel=2,
            )

    if not per_roi_results:
        adata.uns.setdefault("gr", {})["spatial_autocorr"] = {}
        return

    # Align gene index across ROIs
    all_genes = list(dict.fromkeys(g for df in per_roi_results for g in df.index))
    n_roi = len(per_roi_results)
    weights = np.sqrt(np.array(per_roi_n_cells, dtype=float))

    # Combine I (Moran) / C (Geary) via Fisher r-to-z
    stat_col = "I" if mode == "moran" else "C"
    stat_arr = np.full((n_roi, len(all_genes)), np.nan)
    pval_arr = np.full((n_roi, len(all_genes)), 1.0)

    for r, df in enumerate(per_roi_results):
        for g_idx, gene in enumerate(all_genes):
            if gene in df.index:
                stat_arr[r, g_idx] = df.loc[gene, stat_col]
                pval_arr[r, g_idx] = df.loc[gene, "pval_norm"]

    # Fisher r-to-z for Moran I
    clipped = np.clip(stat_arr, -1 + 1e-7, 1 - 1e-7)
    z_transformed = 0.5 * np.log((1 + clipped) / (1 - clipped))
    mean_z = np.nanmean(z_transformed, axis=0)
    stat_combined = np.tanh(mean_z)

    # Stouffer weighted by sqrt(n_cells)
    pval_combined = combine_pvalues(
        pval_arr.reshape(n_roi, 1, len(all_genes)),
        method="stouffer",
        weights=weights,
    )[0]

    pval_fdr = apply_bh_correction(pval_combined)

    result_out = pd.DataFrame(
        {
            stat_col: stat_combined,
            "pval_combined": pval_combined,
            "pval_fdr_bh": pval_fdr,
        },
        index=all_genes,
    )

    adata.uns.setdefault("gr", {})["spatial_autocorr"] = {
        "result": result_out,
        "mode": mode,
    }
