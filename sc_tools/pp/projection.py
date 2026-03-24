"""Subsample and kNN projection for LargeStrategy.

Provides stratified subsampling, kNN index construction (FAISS or BallTree),
and projection of labels, UMAP coordinates, and corrected embeddings from
a representative subsample back to the full dataset.
"""

from __future__ import annotations

import logging
from collections import Counter
from datetime import datetime, timezone

import numpy as np
from anndata import AnnData

logger = logging.getLogger(__name__)

__all__ = [
    "subsample_stratified",
    "build_knn_index",
    "project_labels",
    "project_umap",
    "project_representation",
]


class KNNIndex:
    """Uniform interface wrapping FAISS or sklearn BallTree."""

    def __init__(self, index, backend: str):
        self._index = index
        self.backend = backend

    def query(self, X: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
        """Return (distances, indices) for k nearest neighbors."""
        X = np.ascontiguousarray(X, dtype=np.float32)
        if self.backend == "faiss":
            distances, indices = self._index.search(X, k)
            return distances, indices
        else:  # sklearn BallTree
            distances, indices = self._index.query(X, k=k)
            return distances, indices


def subsample_stratified(
    adata: AnnData,
    n: int = 1_000_000,
    stratify_key: str | None = None,
    random_state: int = 42,
) -> np.ndarray:
    """Proportional stratified subsample with minimum group representation.

    Parameters
    ----------
    adata
        AnnData object to subsample from.
    n
        Target number of cells. If >= n_obs, returns all indices.
    stratify_key
        Column in adata.obs for stratification. If None or not in obs,
        performs simple random subsampling.
    random_state
        Seed for reproducibility.

    Returns
    -------
    np.ndarray of integer indices into adata.obs.
    """
    rng = np.random.RandomState(random_state)  # noqa: NPY002

    if n >= adata.n_obs:
        return np.arange(adata.n_obs)

    # Simple random subsample if no stratification
    if stratify_key is None or stratify_key not in adata.obs.columns:
        return rng.choice(adata.n_obs, size=n, replace=False)

    # Stratified subsample
    groups = adata.obs[stratify_key]
    selected: list[np.ndarray] = []

    # Get value counts only for groups that actually have cells
    value_counts = groups.value_counts()
    non_empty = value_counts[value_counts > 0]

    if len(non_empty) == 0:
        return rng.choice(adata.n_obs, size=n, replace=False)

    total_cells = non_empty.sum()

    # First pass: allocate proportionally with minimum 1 per group
    allocations: dict[str, int] = {}
    for group_name, group_count in non_empty.items():
        proportion = group_count / total_cells
        alloc = max(1, int(round(proportion * n)))
        # Don't allocate more than available
        alloc = min(alloc, group_count)
        allocations[group_name] = alloc

    # Adjust to hit exactly n
    total_alloc = sum(allocations.values())
    if total_alloc > n:
        # Reduce from largest groups
        sorted_groups = sorted(allocations.keys(), key=lambda g: allocations[g], reverse=True)
        for g in sorted_groups:
            if total_alloc <= n:
                break
            reduce = min(allocations[g] - 1, total_alloc - n)
            allocations[g] -= reduce
            total_alloc -= reduce
    elif total_alloc < n:
        # Add to groups that have room
        sorted_groups = sorted(allocations.keys(), key=lambda g: allocations[g], reverse=True)
        for g in sorted_groups:
            if total_alloc >= n:
                break
            room = non_empty[g] - allocations[g]
            add = min(room, n - total_alloc)
            allocations[g] += add
            total_alloc += add

    # Sample from each group
    for group_name, alloc in allocations.items():
        group_idx = np.where(groups == group_name)[0]
        if alloc >= len(group_idx):
            selected.append(group_idx)
        else:
            chosen = rng.choice(group_idx, size=alloc, replace=False)
            selected.append(chosen)

    return np.concatenate(selected)


def build_knn_index(
    adata_sub: AnnData,
    use_rep: str = "X_pca",
) -> KNNIndex:
    """Build kNN index on subsample representation. FAISS > BallTree.

    Parameters
    ----------
    adata_sub
        AnnData with the representation to index.
    use_rep
        Key in obsm for the representation vectors.

    Returns
    -------
    KNNIndex with a uniform query interface.
    """
    data = np.ascontiguousarray(adata_sub.obsm[use_rep], dtype=np.float32)

    # Try FAISS first
    try:
        import faiss  # type: ignore[import-untyped]

        # Try GPU
        try:
            res = faiss.StandardGpuResources()
            index_flat = faiss.IndexFlatL2(data.shape[1])
            index = faiss.index_cpu_to_gpu(res, 0, index_flat)
            logger.info("Built FAISS-GPU kNN index (%d vectors, %d dims)", *data.shape)
        except Exception:
            index = faiss.IndexFlatL2(data.shape[1])
            logger.info("Built FAISS-CPU kNN index (%d vectors, %d dims)", *data.shape)

        index.add(data)
        return KNNIndex(index, backend="faiss")

    except ImportError:
        logger.info("FAISS not available, falling back to sklearn BallTree")

    # Fallback to BallTree
    from sklearn.neighbors import BallTree

    tree = BallTree(data)
    logger.info("Built BallTree kNN index (%d vectors, %d dims)", *data.shape)
    return KNNIndex(tree, backend="sklearn")


def _weighted_average(
    values: np.ndarray,
    distances: np.ndarray,
) -> np.ndarray:
    """Compute distance-weighted average, handling zero distances.

    Parameters
    ----------
    values
        Array of shape (n_queries, k, d) or (n_queries, k).
    distances
        Array of shape (n_queries, k).

    Returns
    -------
    Weighted average of shape (n_queries, d) or (n_queries,).
    """
    # Avoid division by zero: where distance is 0, use a large weight
    eps = 1e-10
    weights = 1.0 / (distances + eps)
    # Normalize weights per query
    weight_sum = weights.sum(axis=1, keepdims=True)
    weights = weights / weight_sum

    if values.ndim == 3:
        # (n, k, d) case
        return np.einsum("nk,nkd->nd", weights, values)
    else:
        # (n, k) case
        return (weights * values).sum(axis=1)


def project_labels(
    adata: AnnData,
    knn_index: KNNIndex,
    adata_sub: AnnData,
    label_key: str = "leiden",
    k: int = 30,
) -> None:
    """kNN majority vote label transfer with confidence.

    Stores ``{label_key}_projected`` and ``{label_key}_confidence`` in
    ``adata.obs``.
    """
    query_data = np.ascontiguousarray(adata.obsm["X_pca"], dtype=np.float32)
    distances, indices = knn_index.query(query_data, k=k)

    sub_labels = adata_sub.obs[label_key].values

    projected = []
    confidence = []
    for i in range(adata.n_obs):
        neighbor_labels = sub_labels[indices[i]]
        counts = Counter(neighbor_labels)
        most_common_label, most_common_count = counts.most_common(1)[0]
        projected.append(most_common_label)
        confidence.append(most_common_count / k)

    adata.obs[f"{label_key}_projected"] = projected
    adata.obs[f"{label_key}_confidence"] = confidence


def project_umap(
    adata: AnnData,
    knn_index: KNNIndex,
    adata_sub: AnnData,
    k: int = 30,
) -> None:
    """Weighted average UMAP coordinate transfer.

    Overwrites ``adata.obsm["X_umap"]``.
    """
    query_data = np.ascontiguousarray(adata.obsm["X_pca"], dtype=np.float32)
    distances, indices = knn_index.query(query_data, k=k)

    sub_umap = adata_sub.obsm["X_umap"]  # (n_sub, 2)

    # Gather neighbor UMAP coords: (n_obs, k, 2)
    neighbor_umap = sub_umap[indices]

    adata.obsm["X_umap"] = _weighted_average(neighbor_umap, distances)


def project_representation(
    adata: AnnData,
    knn_index: KNNIndex,
    adata_sub: AnnData,
    rep_key: str = "X_pca_harmony",
    k: int = 30,
) -> None:
    """Project corrected embedding back to full data.

    Stores the result in ``adata.obsm[rep_key]``.
    """
    query_data = np.ascontiguousarray(adata.obsm["X_pca"], dtype=np.float32)
    distances, indices = knn_index.query(query_data, k=k)

    sub_rep = adata_sub.obsm[rep_key]  # (n_sub, d)

    # Gather neighbor representations: (n_obs, k, d)
    neighbor_rep = sub_rep[indices]

    adata.obsm[rep_key] = _weighted_average(neighbor_rep, distances)


def _store_projection_info(
    adata: AnnData,
    ctx,  # SubsampleContext
    strategy_name: str = "large",
) -> None:
    """Write projection metadata to adata.uns and obs.

    Parameters
    ----------
    adata
        Full AnnData object.
    ctx
        SubsampleContext with subsample configuration.
    strategy_name
        Name of the strategy used.
    """
    adata.uns["_projection_info"] = {
        "strategy": strategy_name,
        "subsample_n": ctx.subsample_n,
        "projection_k": ctx.projection_k,
        "random_state": ctx.random_state,
        "stratify_key": ctx.stratify_key,
        "use_rep": ctx.use_rep,
        "timestamp": datetime.now(tz=timezone.utc).isoformat(),  # noqa: UP017
    }

    # Boolean mask marking representative cells
    mask = np.zeros(adata.n_obs, dtype=bool)
    mask[ctx.subsample_idx] = True
    adata.obs["_sc_tools_is_representative"] = mask
