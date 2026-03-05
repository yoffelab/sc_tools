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
    "run_integration_benchmark",
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
    celltype_key: str | None = None,
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
        Column in ``adata.obs`` with cell type labels. If ``None`` or not
        present in ``adata.obs``, bio conservation metrics are skipped and
        only batch removal metrics are returned.
    use_scib
        ``"auto"`` (default): use scib-metrics if available, else sklearn.
        ``"scib"``: require scib-metrics (error if missing).
        ``"sklearn"``: force sklearn fallbacks.

    Returns
    -------
    Dict with batch removal and bio conservation metric values, all in [0, 1].
    When *celltype_key* is ``None``, bio metrics are omitted.
    """
    if embedding_key not in adata.obsm:
        raise KeyError(f"Embedding {embedding_key!r} not found in adata.obsm")
    if batch_key not in adata.obs.columns:
        raise KeyError(f"Batch key {batch_key!r} not found in adata.obs")

    # Resolve celltype availability
    has_celltype = celltype_key is not None and celltype_key in adata.obs.columns

    X = np.asarray(adata.obsm[embedding_key])
    batch = np.asarray(adata.obs[batch_key])
    celltype = np.asarray(adata.obs[celltype_key]) if has_celltype else None

    should_use_scib = (use_scib == "scib") or (use_scib == "auto" and _HAS_SCIB)
    if use_scib == "scib" and not _HAS_SCIB:
        raise ImportError(
            "scib-metrics is required but not installed. "
            "Install with: pip install sc-tools[benchmark]"
        )

    metrics: dict[str, float] = {}

    if should_use_scib and has_celltype:
        metrics.update(_compute_scib_metrics(adata, embedding_key, batch_key, celltype_key))
    elif should_use_scib and not has_celltype:
        # Batch-only via scib path
        metrics.update(_compute_scib_metrics_batch_only(adata, embedding_key, batch_key))
    else:
        logger.info("scib-metrics not available, using sklearn fallbacks")
        # Batch removal metrics
        metrics["asw_batch"] = _asw_batch_sklearn(X, batch)
        metrics["pcr"] = _pcr_sklearn(X, batch)
        if has_celltype:
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


def _compute_scib_metrics_batch_only(
    adata: AnnData,
    embedding_key: str,
    batch_key: str,
) -> dict[str, float]:
    """Compute batch-only metrics using scib-metrics (no celltype)."""
    import scib_metrics

    X = np.asarray(adata.obsm[embedding_key])
    batch = np.asarray(adata.obs[batch_key])

    metrics: dict[str, float] = {}

    try:
        metrics["asw_batch"] = float(scib_metrics.silhouette_batch(X, batch, np.zeros(len(batch))))
    except Exception:
        metrics["asw_batch"] = _asw_batch_sklearn(X, batch)

    try:
        metrics["pcr"] = float(scib_metrics.pcr_comparison(X, X, batch))
    except Exception:
        metrics["pcr"] = _pcr_sklearn(X, batch)

    try:
        metrics["ilisi"] = float(scib_metrics.ilisi_knn(X, batch))
    except Exception:
        logger.debug("iLISI computation failed, skipping")

    try:
        metrics["kbet"] = float(scib_metrics.kbet(X, batch))
    except Exception:
        logger.debug("kBET computation failed, skipping")

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
    celltype_key: str | None = None,
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
        Column in ``obs`` with cell type labels. If ``None`` or not
        present in ``adata.obs``, bio conservation metrics are skipped.
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

    # TODO: When unintegrated baseline wins, consider emitting a warning.
    # This usually indicates insufficient batch effect or benchmark misconfiguration.

    return df


# ---------------------------------------------------------------------------
# Integration Benchmark Orchestrator
# ---------------------------------------------------------------------------

# Default methods per modality category
_PROTEIN_METHODS = ["harmony", "bbknn", "combat", "scanorama", "cytovi", "pca"]
_TRANSCRIPTOMIC_METHODS = ["harmony", "bbknn", "combat", "scanorama", "scvi", "pca"]

# Modalities that use protein-based integration
_PROTEIN_MODALITIES = {"imc"}


def run_integration_benchmark(
    adata: AnnData,
    modality: str = "visium",
    batch_key: str = "library_id",
    celltype_key: str | None = None,
    methods: list[str] | None = None,
    use_gpu: str | bool = "auto",
    max_epochs: int = 200,
    use_scib: str = "auto",
) -> tuple[AnnData, pd.DataFrame]:
    """Run multiple integration methods and benchmark them.

    Orchestrates integration across methods appropriate for the given
    modality, computes quality metrics for each, and returns the
    AnnData with all embeddings plus a comparison DataFrame.

    Parameters
    ----------
    adata
        AnnData with raw counts (for VAE methods) or normalized data.
        Modified in place with embeddings added to ``obsm``.
    modality
        Data modality (determines default method set and normalization).
    batch_key
        Column in ``adata.obs`` for batch correction.
    celltype_key
        Column in ``adata.obs`` with cell type labels. If provided and
        ``scib-metrics`` is available, the ``Benchmarker`` class is used.
    methods
        List of method names to run. If ``None``, uses modality defaults.
        Valid names: ``harmony``, ``bbknn``, ``combat``, ``scanorama``,
        ``scvi``, ``scanvi``, ``cytovi``, ``pca``.
    use_gpu
        GPU setting for VAE methods.
    max_epochs
        Maximum epochs for VAE methods.
    use_scib
        Passed to metric computation.

    Returns
    -------
    tuple[AnnData, pd.DataFrame]
        The AnnData with all embeddings, and the comparison DataFrame
        sorted by ``overall_score``.
    """
    import scanpy as sc

    is_protein = modality in _PROTEIN_MODALITIES

    if methods is None:
        methods = _PROTEIN_METHODS if is_protein else _TRANSCRIPTOMIC_METHODS

    if batch_key not in adata.obs.columns:
        raise ValueError(f"batch_key '{batch_key}' not in adata.obs")

    # Ensure PCA exists (needed for harmony, bbknn, combat, baseline)
    if "X_pca" not in adata.obsm:
        n_comps = min(50, adata.n_vars - 1, adata.n_obs - 1)
        sc.tl.pca(adata, n_comps=n_comps)
        logger.info("Computed PCA with %d components", n_comps)

    embedding_keys: dict[str, str] = {}

    for method in methods:
        try:
            if method == "harmony":
                from sc_tools.pp.integrate import run_harmony

                run_harmony(adata, batch_key=batch_key)
                embedding_keys["Harmony"] = "X_pca_harmony"

            elif method == "bbknn":
                from sc_tools.pp.integrate import run_bbknn

                run_bbknn(adata, batch_key=batch_key)
                embedding_keys["BBKNN"] = "X_umap_bbknn"

            elif method == "combat":
                from sc_tools.pp.integrate import run_combat

                run_combat(adata, batch_key=batch_key)
                embedding_keys["ComBat"] = "X_pca_combat"

            elif method == "scanorama":
                from sc_tools.pp.integrate import run_scanorama

                run_scanorama(adata, batch_key=batch_key)
                embedding_keys["Scanorama"] = "X_scanorama"

            elif method == "scvi":
                from sc_tools.pp.integrate import run_scvi

                run_scvi(
                    adata,
                    batch_key=batch_key,
                    max_epochs=max_epochs,
                    use_gpu=use_gpu,
                )
                embedding_keys["scVI"] = "X_scVI"

            elif method == "scanvi":
                if celltype_key and celltype_key in adata.obs.columns:
                    from sc_tools.pp.integrate import run_scanvi

                    run_scanvi(
                        adata,
                        batch_key=batch_key,
                        labels_key=celltype_key,
                        max_epochs=max_epochs,
                        use_gpu=use_gpu,
                    )
                    embedding_keys["scANVI"] = "X_scANVI"
                else:
                    logger.info("Skipping scANVI: celltype_key not available")

            elif method == "cytovi":
                from sc_tools.pp.integrate import run_cytovi

                run_cytovi(
                    adata,
                    batch_key=batch_key,
                    max_epochs=max_epochs,
                    use_gpu=use_gpu,
                )
                embedding_keys["CytoVI"] = "X_cytovi"

            elif method == "pca":
                embedding_keys["Unintegrated (PCA)"] = "X_pca"

            else:
                logger.warning("Unknown integration method: %s", method)

        except ImportError as e:
            logger.warning("Skipping %s: %s", method, e)
        except Exception:
            logger.warning("Failed to run %s", method, exc_info=True)

    # Try scib-metrics Benchmarker first
    comparison_df = _try_scib_benchmarker(
        adata,
        embedding_keys,
        batch_key,
        celltype_key,
        use_scib,
    )

    if comparison_df is None:
        # Fall back to our compare_integrations
        comparison_df = compare_integrations(
            adata,
            embedding_keys,
            batch_key,
            celltype_key=celltype_key,
            include_unintegrated=False,
            use_scib=use_scib,
        )

    logger.info("Integration benchmark complete: %d methods evaluated", len(comparison_df))
    return adata, comparison_df


def _try_scib_benchmarker(
    adata: AnnData,
    embedding_keys: dict[str, str],
    batch_key: str,
    celltype_key: str | None,
    use_scib: str,
) -> pd.DataFrame | None:
    """Try to use scib_metrics.benchmark.Benchmarker for comparison.

    Returns None if not available or if celltype_key is missing.
    """
    if use_scib == "sklearn":
        return None
    if celltype_key is None or celltype_key not in adata.obs.columns:
        return None

    try:
        from scib_metrics.benchmark import Benchmarker
    except ImportError:
        return None

    obsm_keys = [v for v in embedding_keys.values() if v in adata.obsm]
    if not obsm_keys:
        return None

    try:
        bm = Benchmarker(
            adata,
            batch_key=batch_key,
            label_key=celltype_key,
            embedding_obsm_keys=obsm_keys,
            n_jobs=-1,
        )
        bm.benchmark()
        df = bm.get_results(min_max_scale=False)

        # Map obsm keys back to method names
        key_to_name = {v: k for k, v in embedding_keys.items()}
        if "Embedding" in df.columns:
            df["method"] = df["Embedding"].map(key_to_name).fillna(df["Embedding"])
        elif df.index.name == "Embedding" or "Embedding" not in df.columns:
            df = df.reset_index()
            if "Embedding" in df.columns:
                df["method"] = df["Embedding"].map(key_to_name).fillna(df["Embedding"])

        # Compute overall score if not present
        if "Total" in df.columns:
            df["overall_score"] = df["Total"]
        elif "overall_score" not in df.columns:
            numeric_cols = df.select_dtypes(include="number").columns
            df["overall_score"] = df[numeric_cols].mean(axis=1)

        df = df.sort_values("overall_score", ascending=False).reset_index(drop=True)
        logger.info("Used scib-metrics Benchmarker for comparison")
        return df

    except Exception:
        logger.debug("scib-metrics Benchmarker failed; falling back", exc_info=True)
        return None
