"""
ScType marker-scoring backend (zero soft dependencies).

Algorithm
---------
1. Load marker DB: ``{celltype: {"positive": [genes], "negative": [genes]}}``.
2. For each cell, compute a score per celltype:
   ``score = sum(expr[pos]) / sqrt(n_pos) - sum(expr[neg]) / sqrt(n_neg)``
   where only genes present in ``adata.var_names`` are used.
3. Assign each cell (or cluster) the celltype with the highest score.
4. Cells (or clusters) below ``min_score`` receive the label ``"Unknown"``.

Reference
---------
Ianevski et al. 2022, Nature Communications.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import scipy.sparse as sp

if TYPE_CHECKING:
    import anndata as ad

from ._base import register_celltype_backend


def _gene_index(genes: list[str], var_upper: dict[str, int]) -> np.ndarray:
    """Return array of var indices for genes present in adata (case-insensitive)."""
    idx = []
    for g in genes:
        i = var_upper.get(g.upper())
        if i is not None:
            idx.append(i)
    return np.array(sorted(set(idx)), dtype=np.intp)


def _compute_scores(
    X: np.ndarray,
    pos_idx: np.ndarray,
    neg_idx: np.ndarray,
) -> np.ndarray:
    """
    Compute per-cell ScType score for one celltype.

    Parameters
    ----------
    X
        Dense (n_obs, n_vars) expression matrix.
    pos_idx
        Column indices for positive markers.
    neg_idx
        Column indices for negative markers.
    """
    n_pos = len(pos_idx)
    n_neg = len(neg_idx)

    pos_score = (
        X[:, pos_idx].sum(axis=1) / np.sqrt(max(n_pos, 1)) if n_pos > 0 else np.zeros(X.shape[0])
    )
    neg_score = (
        X[:, neg_idx].sum(axis=1) / np.sqrt(max(n_neg, 1)) if n_neg > 0 else np.zeros(X.shape[0])
    )

    return pos_score - neg_score


class ScTypeBackend:
    """ScType-style marker scoring; pure numpy, zero soft deps."""

    @staticmethod
    def run(
        adata: ad.AnnData,
        *,
        cluster_key: str = "leiden",
        store_proba: bool = False,
        marker_db: dict,
        use_raw: bool = False,
        scale_factor: float = 1.0,
        min_score: float = 0.0,
        assign_by: str = "cluster",
        **kwargs,
    ) -> tuple[pd.Series, pd.Series, pd.DataFrame | None, dict]:
        """
        Run ScType scoring.

        Parameters
        ----------
        adata
            AnnData with expression in ``X`` (or ``raw.X`` when use_raw=True).
        cluster_key
            Column in obs with cluster labels.
        store_proba
            Whether to return per-celltype score DataFrame as proba.
        marker_db
            Dict mapping celltype name to ``{"positive": [...], "negative": [...]}``.
        use_raw
            Use ``adata.raw.X`` when True.
        scale_factor
            Divide expression by this factor before scoring (e.g. arcsinh cofactor).
        min_score
            Cells (or clusters) with best score below this threshold get "Unknown".
        assign_by
            ``"cluster"`` (default): all cells in a cluster share the cluster-level
            winner.  ``"cell"``: each cell is assigned individually.

        Returns
        -------
        labels, scores, proba_or_None, metadata
        """
        if not marker_db:
            raise ValueError("marker_db is empty; provide at least one celltype entry.")

        # Retrieve expression matrix
        if use_raw and adata.raw is not None:
            X_src = adata.raw.X
            var_names = list(adata.raw.var_names)
        else:
            X_src = adata.X
            var_names = list(adata.var_names)

        # Convert to dense float32
        if sp.issparse(X_src):
            X_dense = X_src.toarray().astype(np.float32)
        else:
            X_dense = np.asarray(X_src, dtype=np.float32)

        if scale_factor != 1.0:
            X_dense = X_dense / scale_factor

        var_upper = {g.upper(): i for i, g in enumerate(var_names)}

        celltypes = list(marker_db.keys())
        n_obs = adata.n_obs
        score_matrix = np.full((n_obs, len(celltypes)), np.nan, dtype=np.float64)

        n_genes_used: dict[str, int] = {}
        for j, ct in enumerate(celltypes):
            entry = marker_db[ct]
            pos_idx = _gene_index(entry.get("positive", []), var_upper)
            neg_idx = _gene_index(entry.get("negative", []), var_upper)
            n_genes_used[ct] = len(pos_idx) + len(neg_idx)
            score_matrix[:, j] = _compute_scores(X_dense, pos_idx, neg_idx)

        score_df = pd.DataFrame(score_matrix, index=adata.obs_names, columns=celltypes)

        if assign_by == "cluster":
            labels_arr, scores_arr = _assign_by_cluster(
                score_df, adata.obs[cluster_key].astype(str), celltypes, min_score
            )
        else:
            labels_arr, scores_arr = _assign_by_cell(score_df, celltypes, min_score)

        labels = pd.Series(labels_arr, index=adata.obs_names, dtype=object)
        scores = pd.Series(scores_arr, index=adata.obs_names, dtype=np.float64)
        proba = score_df if store_proba else None

        meta = {
            "n_celltypes": len(celltypes),
            "assign_by": assign_by,
            "min_score": min_score,
            "n_genes_used": n_genes_used,
        }
        return labels, scores, proba, meta


def _assign_by_cell(
    score_df: pd.DataFrame,
    celltypes: list[str],
    min_score: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Assign each cell individually to the highest-scoring celltype."""
    scores_arr = score_df.values  # (n_obs, n_celltypes)
    best_idx = np.nanargmax(scores_arr, axis=1)
    best_scores = scores_arr[np.arange(len(best_idx)), best_idx]

    labels_out = np.array([celltypes[i] for i in best_idx], dtype=object)
    labels_out[best_scores < min_score] = "Unknown"

    return labels_out, np.where(best_scores < min_score, 0.0, best_scores)


def _assign_by_cluster(
    score_df: pd.DataFrame,
    cluster_series: pd.Series,
    celltypes: list[str],
    min_score: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Assign each cluster to its highest-scoring celltype; broadcast to cells."""
    clusters = cluster_series.values
    unique_clusters = np.unique(clusters)

    cell_labels = np.full(len(score_df), "Unknown", dtype=object)
    cell_scores = np.zeros(len(score_df), dtype=np.float64)

    for cl in unique_clusters:
        mask = clusters == cl
        cluster_scores = score_df.values[mask].mean(axis=0)  # mean score per celltype
        best_j = int(np.nanargmax(cluster_scores))
        best_score = float(cluster_scores[best_j])
        if best_score >= min_score:
            cell_labels[mask] = celltypes[best_j]
            cell_scores[mask] = best_score

    return cell_labels, cell_scores


register_celltype_backend("sctype", ScTypeBackend)
