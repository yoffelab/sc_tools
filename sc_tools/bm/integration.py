"""
Batch correction quality metrics (scib-style).

Computes batch-removal and bio-conservation metrics for integration
method comparison. Uses scib-metrics when available, falls back to
sklearn for core metrics (ASW, ARI, NMI, PCR).
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from anndata import AnnData

__all__ = [
    "compute_integration_metrics",
    "compute_composite_score",
    "compare_integrations",
]

logger = logging.getLogger(__name__)

# Check for scib-metrics availability
try:
    import scib_metrics  # noqa: F401

    _HAS_SCIB = True
except ImportError:
    _HAS_SCIB = False


# ---------------------------------------------------------------------------
# Metric computation helpers (sklearn fallbacks)
# ---------------------------------------------------------------------------


def _asw_batch_sklearn(X: np.ndarray, batch: np.ndarray) -> float:
    """ASW batch removal: 1 - |silhouette_score| (higher = better mixing)."""
    from sklearn.metrics import silhouette_score

    n_unique = len(np.unique(batch))
    if n_unique < 2 or len(X) < n_unique + 1:
        return 0.0
    score = silhouette_score(X, batch, sample_size=min(len(X), 5000), random_state=42)
    return float(1 - abs(score))


def _asw_celltype_sklearn(X: np.ndarray, celltype: np.ndarray) -> float:
    """ASW bio conservation: silhouette_score on cell types (higher = better separation)."""
    from sklearn.metrics import silhouette_score

    n_unique = len(np.unique(celltype))
    if n_unique < 2 or len(X) < n_unique + 1:
        return 0.0
    score = silhouette_score(X, celltype, sample_size=min(len(X), 5000), random_state=42)
    # Scale from [-1, 1] to [0, 1]
    return float((score + 1) / 2)


def _pcr_sklearn(X: np.ndarray, batch: np.ndarray) -> float:
    """Principal component regression: 1 - R^2 of batch on PCs."""
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LinearRegression
    from sklearn.preprocessing import LabelEncoder

    if len(np.unique(batch)) < 2:
        return 0.0

    n_components = min(50, X.shape[1], X.shape[0] - 1)
    if n_components < 1:
        return 0.0

    pca = PCA(n_components=n_components)
    pcs = pca.fit_transform(X)

    le = LabelEncoder()
    batch_encoded = le.fit_transform(batch).reshape(-1, 1)

    lr = LinearRegression()
    lr.fit(batch_encoded, pcs)
    r2 = lr.score(batch_encoded, pcs)
    return float(1 - max(0, r2))


def _ari_sklearn(X: np.ndarray, celltype: np.ndarray, resolution: float = 1.0) -> float:
    """ARI: cluster (Leiden) vs cell type labels."""
    from sklearn.metrics import adjusted_rand_score

    clusters = _leiden_cluster(X, resolution=resolution)
    return float(adjusted_rand_score(celltype, clusters))


def _nmi_sklearn(X: np.ndarray, celltype: np.ndarray, resolution: float = 1.0) -> float:
    """NMI: cluster (Leiden) vs cell type labels."""
    from sklearn.metrics import normalized_mutual_info_score

    clusters = _leiden_cluster(X, resolution=resolution)
    return float(normalized_mutual_info_score(celltype, clusters, average_method="arithmetic"))


def _graph_connectivity(X: np.ndarray, celltype: np.ndarray) -> float:
    """Graph connectivity: fraction of cells in largest connected component per cell type."""
    import scanpy as sc

    tmp = AnnData(X=np.zeros((len(celltype), 1)))
    tmp.obsm["X_emb"] = X
    sc.pp.neighbors(tmp, use_rep="X_emb", n_neighbors=15)

    from scipy.sparse.csgraph import connected_components

    adj = tmp.obsp["connectivities"]
    unique_types = np.unique(celltype)
    fractions = []
    for ct in unique_types:
        ct_mask = celltype == ct
        ct_indices = np.where(ct_mask)[0]
        if len(ct_indices) < 2:
            fractions.append(1.0)
            continue
        sub_adj = adj[ct_indices][:, ct_indices]
        n_components, labels = connected_components(sub_adj, directed=False)
        largest_cc = np.bincount(labels).max()
        fractions.append(largest_cc / len(ct_indices))

    return float(np.mean(fractions))


def _leiden_cluster(X: np.ndarray, resolution: float = 1.0) -> np.ndarray:
    """Run Leiden clustering on an embedding."""
    import scanpy as sc

    tmp = AnnData(X=np.zeros((X.shape[0], 1)))
    tmp.obsm["X_emb"] = X
    sc.pp.neighbors(tmp, use_rep="X_emb", n_neighbors=15)
    sc.tl.leiden(tmp, resolution=resolution, key_added="leiden")
    return tmp.obs["leiden"].values


# ---------------------------------------------------------------------------
# Main API
# ---------------------------------------------------------------------------


def compute_integration_metrics(
    adata: AnnData,
    embedding_key: str,
    batch_key: str,
    celltype_key: str,
    use_scib: str = "auto",
) -> dict[str, float]:
    """Compute integration quality metrics for one embedding.

    Parameters
    ----------
    adata
        AnnData object with embeddings in ``obsm``.
    embedding_key
        Key in ``adata.obsm`` (e.g. ``"X_scVI"``).
    batch_key
        Column in ``adata.obs`` with batch labels.
    celltype_key
        Column in ``adata.obs`` with cell type labels.
    use_scib
        ``"auto"`` (default): use scib-metrics if available, else sklearn.
        ``"scib"``: require scib-metrics (error if missing).
        ``"sklearn"``: force sklearn fallbacks.

    Returns
    -------
    Dict with batch removal and bio conservation metric values, all in [0, 1].
    """
    if embedding_key not in adata.obsm:
        raise KeyError(f"Embedding {embedding_key!r} not found in adata.obsm")
    if batch_key not in adata.obs.columns:
        raise KeyError(f"Batch key {batch_key!r} not found in adata.obs")
    if celltype_key not in adata.obs.columns:
        raise KeyError(f"Cell type key {celltype_key!r} not found in adata.obs")

    X = np.asarray(adata.obsm[embedding_key])
    batch = np.asarray(adata.obs[batch_key])
    celltype = np.asarray(adata.obs[celltype_key])

    should_use_scib = (use_scib == "scib") or (use_scib == "auto" and _HAS_SCIB)
    if use_scib == "scib" and not _HAS_SCIB:
        raise ImportError(
            "scib-metrics is required but not installed. "
            "Install with: pip install sc-tools[benchmark]"
        )

    metrics: dict[str, float] = {}

    if should_use_scib:
        metrics.update(_compute_scib_metrics(adata, embedding_key, batch_key, celltype_key))
    else:
        logger.info("scib-metrics not available, using sklearn fallbacks")
        # Batch removal metrics
        metrics["asw_batch"] = _asw_batch_sklearn(X, batch)
        metrics["pcr"] = _pcr_sklearn(X, batch)
        metrics["graph_connectivity"] = _graph_connectivity(X, celltype)

        # Bio conservation metrics
        metrics["asw_celltype"] = _asw_celltype_sklearn(X, celltype)
        metrics["ari"] = _ari_sklearn(X, celltype)
        metrics["nmi"] = _nmi_sklearn(X, celltype)

    # Clamp to [0, 1]
    for k, v in metrics.items():
        metrics[k] = float(np.clip(v, 0.0, 1.0))

    return metrics


def _compute_scib_metrics(
    adata: AnnData,
    embedding_key: str,
    batch_key: str,
    celltype_key: str,
) -> dict[str, float]:
    """Compute metrics using scib-metrics library."""
    import scib_metrics

    X = np.asarray(adata.obsm[embedding_key])
    batch = np.asarray(adata.obs[batch_key])
    celltype = np.asarray(adata.obs[celltype_key])

    metrics: dict[str, float] = {}

    # Batch removal
    try:
        metrics["asw_batch"] = float(scib_metrics.silhouette_batch(X, batch, celltype))
    except Exception:
        metrics["asw_batch"] = _asw_batch_sklearn(X, batch)

    try:
        metrics["pcr"] = float(scib_metrics.pcr_comparison(X, X, batch))
    except Exception:
        metrics["pcr"] = _pcr_sklearn(X, batch)

    metrics["graph_connectivity"] = _graph_connectivity(X, celltype)

    try:
        metrics["ilisi"] = float(scib_metrics.ilisi_knn(X, batch))
    except Exception:
        logger.debug("iLISI computation failed, skipping")

    try:
        metrics["kbet"] = float(scib_metrics.kbet(X, batch))
    except Exception:
        logger.debug("kBET computation failed, skipping")

    # Bio conservation
    try:
        metrics["asw_celltype"] = float(scib_metrics.silhouette_label(X, celltype))
    except Exception:
        metrics["asw_celltype"] = _asw_celltype_sklearn(X, celltype)

    metrics["ari"] = _ari_sklearn(X, celltype)
    metrics["nmi"] = _nmi_sklearn(X, celltype)

    try:
        metrics["clisi"] = float(scib_metrics.clisi_knn(X, celltype))
    except Exception:
        logger.debug("cLISI computation failed, skipping")

    try:
        metrics["isolated_label_f1"] = float(scib_metrics.isolated_labels(X, celltype, batch))
    except Exception:
        logger.debug("Isolated label F1 computation failed, skipping")

    return metrics


def compute_composite_score(
    metrics: dict[str, float],
    batch_weight: float = 0.4,
    bio_weight: float = 0.6,
) -> dict[str, float]:
    """Compute composite integration score from individual metrics.

    Parameters
    ----------
    metrics
        Output from ``compute_integration_metrics``.
    batch_weight
        Weight for batch removal score (default 0.4).
    bio_weight
        Weight for bio conservation score (default 0.6).

    Returns
    -------
    Dict with ``batch_score``, ``bio_score``, ``overall_score``.
    """
    batch_keys = ["asw_batch", "pcr", "graph_connectivity", "ilisi", "kbet"]
    bio_keys = ["asw_celltype", "ari", "nmi", "clisi", "isolated_label_f1"]

    batch_vals = [metrics[k] for k in batch_keys if k in metrics]
    bio_vals = [metrics[k] for k in bio_keys if k in metrics]

    batch_score = float(np.mean(batch_vals)) if batch_vals else 0.0
    bio_score = float(np.mean(bio_vals)) if bio_vals else 0.0

    overall = batch_weight * batch_score + bio_weight * bio_score

    return {
        "batch_score": batch_score,
        "bio_score": bio_score,
        "overall_score": float(overall),
    }


def compare_integrations(
    adata: AnnData,
    embeddings: dict[str, str],
    batch_key: str,
    celltype_key: str,
    batch_weight: float = 0.4,
    bio_weight: float = 0.6,
    include_unintegrated: bool = True,
    use_scib: str = "auto",
) -> pd.DataFrame:
    """Compare multiple integration methods side-by-side.

    Parameters
    ----------
    adata
        AnnData with multiple embeddings in ``obsm``.
    embeddings
        Dict mapping method name to ``obsm`` key
        (e.g. ``{"scVI": "X_scVI", "Harmony": "X_pca_harmony"}``).
    batch_key
        Column in ``obs`` with batch labels.
    celltype_key
        Column in ``obs`` with cell type labels.
    batch_weight
        Weight for batch removal in composite score.
    bio_weight
        Weight for bio conservation in composite score.
    include_unintegrated
        If True and ``"X_pca"`` exists, add unintegrated PCA as baseline.
    use_scib
        Passed to ``compute_integration_metrics``.

    Returns
    -------
    DataFrame with rows=methods, sorted by ``overall_score`` descending.
    """
    if include_unintegrated and "X_pca" in adata.obsm and "Unintegrated" not in embeddings:
        embeddings = {"Unintegrated": "X_pca", **embeddings}

    rows = []
    for name, key in embeddings.items():
        if key not in adata.obsm:
            logger.warning("Embedding %r (%s) not in adata.obsm, skipping", name, key)
            continue

        metrics = compute_integration_metrics(
            adata, key, batch_key, celltype_key, use_scib=use_scib
        )
        composite = compute_composite_score(metrics, batch_weight, bio_weight)

        row = {"method": name, "embedding_key": key}
        row.update(metrics)
        row.update(composite)
        rows.append(row)

    df = pd.DataFrame(rows)
    if len(df) > 0:
        df = df.sort_values("overall_score", ascending=False).reset_index(drop=True)

    return df
