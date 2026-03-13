"""
Matrix unification and p-value combination helpers for sc_tools.gr.
"""

from __future__ import annotations

import numpy as np
import scipy.stats as ss
from statsmodels.stats.multitest import multipletests


def unify_matrices(
    per_roi_mats: list[np.ndarray],
    per_roi_cats: list[list[str]],
    fill_value: float = np.nan,
) -> tuple[np.ndarray, list[str]]:
    """
    Unify a list of per-ROI 2D square matrices into a (n_roi, N_global, N_global) array.

    Categories present in some ROIs but absent in others are filled with fill_value.

    Parameters
    ----------
    per_roi_mats
        List of 2D arrays, one per ROI.  Shape: (n_cats_roi, n_cats_roi).
    per_roi_cats
        List of category lists (same order as matrix axes) for each ROI.
    fill_value
        Value to use for absent categories.

    Returns
    -------
    unified
        Array of shape (n_roi, N_global, N_global).
    global_cats
        Ordered list of all categories across all ROIs.
    """
    # Build ordered global category list (preserving first-seen order)
    seen: dict[str, None] = {}
    for cats in per_roi_cats:
        for c in cats:
            seen[c] = None
    global_cats: list[str] = list(seen.keys())
    n_global = len(global_cats)
    cat_to_idx = {c: i for i, c in enumerate(global_cats)}

    n_roi = len(per_roi_mats)
    unified = np.full((n_roi, n_global, n_global), fill_value, dtype=float)

    for r, (mat, cats) in enumerate(zip(per_roi_mats, per_roi_cats, strict=False)):
        for i, ci in enumerate(cats):
            gi = cat_to_idx[ci]
            for j, cj in enumerate(cats):
                gj = cat_to_idx[cj]
                unified[r, gi, gj] = mat[i, j]

    return unified, global_cats


def unify_dataframes(
    per_roi_dfs: list,
    fill_value: float = np.nan,
) -> list:
    """
    Align a list of per-ROI DataFrames to a common index (union of all indices).

    Parameters
    ----------
    per_roi_dfs
        List of DataFrames with the same columns but possibly different row indices
        (one index entry per cell type).
    fill_value
        Value for missing rows.

    Returns
    -------
    List of aligned DataFrames, all with the same index in global order.
    """
    # Union of indices, preserving first-seen order
    seen: dict[str, None] = {}
    for df in per_roi_dfs:
        for idx in df.index:
            seen[idx] = None
    global_index = list(seen.keys())

    aligned = []
    for df in per_roi_dfs:
        reindexed = df.reindex(global_index, fill_value=fill_value)
        aligned.append(reindexed)
    return aligned


def combine_pvalues(
    pval_array: np.ndarray,
    method: str = "stouffer",
    weights: np.ndarray | None = None,
) -> np.ndarray:
    """
    Combine p-values across ROIs (axis 0) using Stouffer or Fisher method.

    Parameters
    ----------
    pval_array
        Array of shape (n_roi, ...) with p-values along axis 0.
    method
        "stouffer" or "fisher".
    weights
        Optional weight array broadcastable to (n_roi, ...).
        Only used for Stouffer.  Defaults to uniform weights.

    Returns
    -------
    combined
        Array of shape (...) with combined p-values.
    """
    arr = np.asarray(pval_array, dtype=float)
    # Clip to avoid log(0) / inf in transformations
    arr = np.clip(arr, 1e-300, 1 - 1e-7)

    n_roi = arr.shape[0]
    rest_shape = arr.shape[1:]

    if method == "stouffer":
        z_scores = ss.norm.ppf(1.0 - arr)  # inverse survival: Phi^{-1}(1-p)
        if weights is None:
            w = np.ones(n_roi)
        else:
            w = np.asarray(weights, dtype=float)
        # Broadcast weights to full array shape
        w_broad = w.reshape((n_roi,) + (1,) * len(rest_shape))
        z_combined = np.nansum(w_broad * z_scores, axis=0) / np.sqrt(
            np.nansum(w_broad**2 * np.isfinite(z_scores), axis=0)
        )
        combined = ss.norm.sf(z_combined)

    elif method == "fisher":
        chi2 = -2.0 * np.nansum(np.log(arr), axis=0)
        k = np.sum(np.isfinite(arr), axis=0)
        combined = ss.chi2.sf(chi2, df=2 * k)

    else:
        raise ValueError(f"Unknown method '{method}'. Choose 'stouffer' or 'fisher'.")

    return np.clip(combined, 0.0, 1.0)


def cross_roi_zscore(unified_array: np.ndarray) -> np.ndarray:
    """
    Compute z-scores across ROIs (axis 0) using nanmean / nanstd.

    Requires at least 3 non-NaN ROIs per element; otherwise returns NaN.

    Parameters
    ----------
    unified_array
        Array of shape (n_roi, ...).

    Returns
    -------
    z
        Array of same shape as unified_array with z-scores.
    """
    arr = np.asarray(unified_array, dtype=float)
    valid_count = np.sum(np.isfinite(arr), axis=0, keepdims=True)  # (1, ...)
    mean = np.nanmean(arr, axis=0, keepdims=True)
    std = np.nanstd(arr, axis=0, keepdims=True)

    z = np.where(valid_count >= 3, (arr - mean) / np.where(std == 0, np.nan, std), np.nan)
    return z


def apply_bh_correction(pval_matrix: np.ndarray) -> np.ndarray:
    """
    Apply Benjamini-Hochberg FDR correction to all entries of a flat or 2D p-value array.

    Parameters
    ----------
    pval_matrix
        Array of any shape containing p-values.

    Returns
    -------
    fdr
        Array of same shape with BH-corrected q-values.
    """
    orig_shape = pval_matrix.shape
    flat = pval_matrix.ravel()
    finite_mask = np.isfinite(flat)
    fdr_flat = np.full_like(flat, np.nan)
    if finite_mask.any():
        _, q, _, _ = multipletests(flat[finite_mask], method="fdr_bh")
        fdr_flat[finite_mask] = q
    return fdr_flat.reshape(orig_shape)
