"""
Batch correction quality metrics (scib-style).

Computes batch-removal and bio-conservation metrics for integration
method comparison. Uses scib-metrics when available, falls back to
sklearn for core metrics (ASW, ARI, NMI, PCR).
"""

from __future__ import annotations

import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData

__all__ = [
    "compute_integration_metrics",
    "compute_composite_score",
    "compare_integrations",
    "run_integration_benchmark",
    "run_full_integration_workflow",
    "_mask_nan_rows",
    "_load_embedding_h5py",
]

logger = logging.getLogger(__name__)

# Check for scib-metrics availability
try:
    import scib_metrics  # noqa: F401

    _HAS_SCIB = True
except ImportError:
    _HAS_SCIB = False


# ---------------------------------------------------------------------------
# NaN masking helper
# ---------------------------------------------------------------------------


def _mask_nan_rows(X: np.ndarray, *arrays: np.ndarray) -> tuple[np.ndarray, ...]:
    """Remove rows where X has any NaN; apply same mask to companion arrays."""
    valid = ~np.isnan(X).any(axis=1)
    n_dropped = int((~valid).sum())
    if n_dropped > 0:
        logger.warning(
            "Dropped %d cells with NaN embeddings (%d remain)",
            n_dropped,
            int(valid.sum()),
        )
    return (X[valid], *(a[valid] for a in arrays))


# ---------------------------------------------------------------------------
# h5py embedding loader
# ---------------------------------------------------------------------------


def _load_embedding_h5py(
    path: str | Path,
    obsm_key: str,
    obs_keys: list[str],
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Load embedding + obs columns from h5ad via h5py (no full AnnData).

    Parameters
    ----------
    path
        Path to an h5ad file.
    obsm_key
        Key in ``obsm`` group (e.g. ``"X_scVI"``).
    obs_keys
        List of obs column names to load (e.g. ``["batch", "celltype"]``).

    Returns
    -------
    tuple of (embedding_array, obs_dict)
        embedding_array has shape (n_obs, n_latent).
        obs_dict maps column name to 1-D array.
    """
    import h5py

    with h5py.File(path, "r") as f:
        # Load embedding
        obsm_path = f"obsm/{obsm_key}"
        if obsm_path not in f:
            raise KeyError(f"Embedding {obsm_key!r} not found in {path}")
        embedding = f[obsm_path][:]

        # Load obs columns
        obs_data: dict[str, np.ndarray] = {}
        for key in obs_keys:
            obs_path = f"obs/{key}"
            if obs_path not in f:
                raise KeyError(f"obs column {key!r} not found in {path}")
            grp = f[obs_path]
            if isinstance(grp, h5py.Group) and "categories" in grp:
                # Categorical column: reconstruct from codes + categories
                codes = grp["codes"][:]
                cats = grp["categories"][:]
                if cats.dtype.kind in ("O", "S", "U"):
                    cats = cats.astype(str)
                obs_data[key] = cats[codes]
            else:
                # Plain array (numeric or string dataset)
                obs_data[key] = grp[:]

    return embedding, obs_data


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


def _ari_sklearn(
    X: np.ndarray, celltype: np.ndarray, resolution: float = 1.0, random_state: int = 0,
) -> float:
    """ARI: cluster (Leiden) vs cell type labels."""
    from sklearn.metrics import adjusted_rand_score

    clusters = _leiden_cluster(X, resolution=resolution, random_state=random_state)
    return float(adjusted_rand_score(celltype, clusters))


def _nmi_sklearn(
    X: np.ndarray, celltype: np.ndarray, resolution: float = 1.0, random_state: int = 0,
) -> float:
    """NMI: cluster (Leiden) vs cell type labels."""
    from sklearn.metrics import normalized_mutual_info_score

    clusters = _leiden_cluster(X, resolution=resolution, random_state=random_state)
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


def _leiden_cluster(X: np.ndarray, resolution: float = 1.0, random_state: int = 0) -> np.ndarray:
    """Run Leiden clustering on an embedding.

    Parameters
    ----------
    X
        Embedding array (n_cells, n_dims).
    resolution
        Leiden clustering resolution.
    random_state
        Random state for reproducibility (D-14, PRV-05).
    """
    import scanpy as sc

    tmp = AnnData(X=np.zeros((X.shape[0], 1)))
    tmp.obsm["X_emb"] = X
    sc.pp.neighbors(tmp, use_rep="X_emb", n_neighbors=15)
    sc.tl.leiden(tmp, resolution=resolution, key_added="leiden", random_state=random_state)
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
    resolution: float = 1.0,
    random_state: int = 0,
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
        metrics.update(_compute_scib_metrics(adata, embedding_key, batch_key, celltype_key, resolution=resolution, random_state=random_state))
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
            metrics["ari"] = _ari_sklearn(X, celltype, resolution=resolution, random_state=random_state)
            metrics["nmi"] = _nmi_sklearn(X, celltype, resolution=resolution, random_state=random_state)

    # Clamp to [0, 1]
    for k, v in metrics.items():
        metrics[k] = float(np.clip(v, 0.0, 1.0))

    return metrics


def _compute_scib_metrics(
    adata: AnnData,
    embedding_key: str,
    batch_key: str,
    celltype_key: str,
    resolution: float = 1.0,
    random_state: int = 0,
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

    metrics["ari"] = _ari_sklearn(X, celltype, resolution=resolution, random_state=random_state)
    metrics["nmi"] = _nmi_sklearn(X, celltype, resolution=resolution, random_state=random_state)

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
    batch_weight: float = 0.2,
    bio_weight: float = 0.8,
) -> dict[str, float]:
    """Compute composite integration score from individual metrics.

    Parameters
    ----------
    metrics
        Output from ``compute_integration_metrics``.
    batch_weight
        Weight for batch removal score (default 0.2).
    bio_weight
        Weight for bio conservation score (default 0.8).

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
    adata: AnnData | None = None,
    embeddings: dict[str, str] | None = None,
    batch_key: str = "batch",
    celltype_key: str | None = None,
    bio_key: str | None = None,
    batch_weight: float = 0.2,
    bio_weight: float = 0.8,
    include_unintegrated: bool = True,
    use_scib: str = "auto",
    subsample_n: int | None = None,
    seed: int = 42,
    resolution: float = 1.0,
    random_state: int = 0,
    embedding_files: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Compare multiple integration methods side-by-side.

    Parameters
    ----------
    adata
        AnnData with multiple embeddings in ``obsm``. Can be ``None`` when
        ``embedding_files`` provides all methods.
    embeddings
        Dict mapping method name to ``obsm`` key
        (e.g. ``{"scVI": "X_scVI", "Harmony": "X_pca_harmony"}``).
    batch_key
        Column in ``obs`` with batch labels.
    celltype_key
        Column in ``obs`` with cell type labels. If ``None`` or not
        present in ``adata.obs``, bio conservation metrics are skipped.
    bio_key
        Column to use for bio conservation metrics. Defaults to
        ``celltype_key`` when not provided. Allows using any clinically
        relevant variable (e.g. ``"condition"``, ``"disease_status"``).
    batch_weight
        Weight for batch removal in composite score.
    bio_weight
        Weight for bio conservation in composite score.
    include_unintegrated
        If True and ``"X_pca"`` exists, add unintegrated PCA as baseline.
    use_scib
        Passed to ``compute_integration_metrics``.
    subsample_n
        If set, subsample to this many cells (stratified by ``batch_key``)
        before computing metrics.
    seed
        Random seed for subsampling reproducibility.
    resolution
        Leiden clustering resolution for ARI/NMI bio conservation metrics.
    embedding_files
        Dict mapping method name to h5ad file path. Embeddings are loaded
        via h5py without reading the full AnnData, enabling benchmarking
        of datasets too large to fit in memory.

    Returns
    -------
    DataFrame with rows=methods, sorted by ``overall_score`` descending.
    """
    # Validate inputs
    if adata is None and not embedding_files:
        raise ValueError("Either adata or embedding_files must be provided")

    if embeddings is None:
        embeddings = {}

    # Resolve bio_key -- defaults to celltype_key for backwards compatibility
    effective_bio_key = bio_key if bio_key is not None else celltype_key

    if adata is not None and include_unintegrated and "X_pca" in adata.obsm and "Unintegrated" not in embeddings:
        embeddings = {"Unintegrated": "X_pca", **embeddings}

    # Optional subsampling before metric computation
    if adata is not None and subsample_n is not None and subsample_n < adata.n_obs:
        adata = _stratified_subsample(adata, batch_key, n=subsample_n, seed=seed)

    rows = []

    # --- adata-based embeddings ---
    if adata is not None:
        for name, key in embeddings.items():
            if key not in adata.obsm:
                logger.warning("Embedding %r (%s) not in adata.obsm, skipping", name, key)
                continue

            # Per-embedding NaN masking (BM-02)
            X_emb = np.asarray(adata.obsm[key])
            valid = ~np.isnan(X_emb).any(axis=1)
            n_dropped = int((~valid).sum())
            if n_dropped > 0:
                logger.warning("Method %s: dropped %d cells with NaN embeddings", name, n_dropped)
                adata_clean = adata[valid].copy()
            else:
                adata_clean = adata

            metrics = compute_integration_metrics(
                adata_clean, key, batch_key, effective_bio_key, use_scib=use_scib, resolution=resolution, random_state=random_state
            )
            composite = compute_composite_score(metrics, batch_weight, bio_weight)

            row = {"method": name, "embedding_key": key}
            row.update(metrics)
            row.update(composite)
            rows.append(row)

    # --- file-based embeddings (BM-01) ---
    if embedding_files:
        import h5py as _h5py

        obs_keys_needed = [batch_key]
        if effective_bio_key:
            obs_keys_needed.append(effective_bio_key)

        for name, fpath in embedding_files.items():
            try:
                # Discover obsm key from file (use the first/only one)
                with _h5py.File(fpath, "r") as f:
                    available_obsm = list(f["obsm"].keys())
                if not available_obsm:
                    logger.warning("No obsm keys in %s, skipping %s", fpath, name)
                    continue
                obsm_key = available_obsm[0]

                emb_array, obs_arrays = _load_embedding_h5py(fpath, obsm_key, obs_keys_needed)

                # Build minimal AnnData for metric computation
                file_adata = AnnData(X=np.zeros((emb_array.shape[0], 1), dtype=np.float32))
                file_adata.obsm[obsm_key] = emb_array
                for k, v in obs_arrays.items():
                    if v.dtype.kind in ("U", "O", "S"):
                        file_adata.obs[k] = pd.Categorical(v)
                    else:
                        file_adata.obs[k] = v

                # Per-embedding NaN masking
                valid = ~np.isnan(emb_array).any(axis=1)
                n_dropped = int((~valid).sum())
                if n_dropped > 0:
                    logger.warning("Method %s: dropped %d cells with NaN embeddings", name, n_dropped)
                    file_adata = file_adata[valid].copy()

                metrics = compute_integration_metrics(
                    file_adata, obsm_key, batch_key, effective_bio_key, use_scib=use_scib, resolution=resolution, random_state=random_state
                )
                composite = compute_composite_score(metrics, batch_weight, bio_weight)

                row = {"method": name, "embedding_key": obsm_key}
                row.update(metrics)
                row.update(composite)
                rows.append(row)

            except Exception:
                logger.warning("Failed to load embedding from %s for %s", fpath, name, exc_info=True)

    df = pd.DataFrame(rows)
    if len(df) > 0:
        df = df.sort_values("overall_score", ascending=False).reset_index(drop=True)

    # Annotate whether sklearn fallback was actually used so that report consumers
    # can surface a warning.  Fallback is active when scib-metrics is absent OR
    # when the caller explicitly forced use_scib="sklearn".
    _sklearn_forced = use_scib == "sklearn"
    df.attrs["scib_fallback"] = _sklearn_forced or not _HAS_SCIB

    # Store benchmark parameters for provenance (BM-06, D-15)
    df.attrs["benchmark_params"] = {
        "batch_weight": batch_weight,
        "bio_weight": bio_weight,
        "subsample_n": subsample_n,
        "seed": seed,
        "resolution": resolution,
        "random_state": random_state,
        "use_scib": use_scib,
    }

    return df


# ---------------------------------------------------------------------------
# Integration Benchmark Orchestrator
# ---------------------------------------------------------------------------

# Default methods per modality category
_PROTEIN_METHODS = ["harmony", "bbknn", "combat", "scanorama", "cytovi", "pca"]
_TRANSCRIPTOMIC_METHODS = ["harmony", "bbknn", "combat", "scanorama", "scvi", "resolvi", "pca"]

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
    runtimes: dict[str, float] = {}

    for method in methods:
        t0 = time.perf_counter()
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

            elif method == "resolvi":
                from sc_tools.pp.integrate import run_resolvi

                run_resolvi(
                    adata,
                    batch_key=batch_key,
                    max_epochs=max_epochs,
                    use_gpu=use_gpu,
                )
                embedding_keys["resolVI"] = "X_resolvi"

            elif method == "pca":
                embedding_keys["Unintegrated (PCA)"] = "X_pca"

            else:
                logger.warning("Unknown integration method: %s", method)

        except ImportError as e:
            logger.warning("Skipping %s: %s", method, e)
        except Exception:
            logger.warning("Failed to run %s", method, exc_info=True)
        finally:
            runtime = time.perf_counter() - t0
            # Record runtime for any display names added during this iteration
            for display_name in list(embedding_keys):
                if display_name not in runtimes:
                    runtimes[display_name] = runtime

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

    # Add runtime column (BM-05)
    if "method" in comparison_df.columns:
        comparison_df["runtime_s"] = comparison_df["method"].map(runtimes).fillna(0.0)

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


# ---------------------------------------------------------------------------
# Full Integration Workflow
# ---------------------------------------------------------------------------

# Maps method short name -> (integration function name, obsm key)
_METHOD_APPLY: dict[str, tuple[str | None, str]] = {
    "harmony": ("run_harmony", "X_pca_harmony"),
    "bbknn": ("run_bbknn", "X_umap_bbknn"),
    "combat": ("run_combat", "X_pca_combat"),
    "scanorama": ("run_scanorama", "X_scanorama"),
    "scvi": ("run_scvi", "X_scVI"),
    "scanvi": ("run_scanvi", "X_scANVI"),
    "cytovi": ("run_cytovi", "X_cytovi"),
    "resolvi": ("run_resolvi", "X_resolvi"),
    "pca": (None, "X_pca"),
}

# Maps display names (from benchmark results) back to method short names
_DISPLAY_TO_METHOD: dict[str, str] = {
    "Harmony": "harmony",
    "BBKNN": "bbknn",
    "ComBat": "combat",
    "Scanorama": "scanorama",
    "scVI": "scvi",
    "scANVI": "scanvi",
    "CytoVI": "cytovi",
    "resolVI": "resolvi",
    "Unintegrated (PCA)": "pca",
    "Unintegrated": "pca",
}


def _stratified_subsample(
    adata: AnnData,
    key: str,
    n: int | None = None,
    fraction: float | None = None,
    seed: int = 42,
) -> AnnData:
    """Subsample adata stratified by key, preserving proportions."""
    rng = np.random.RandomState(seed)

    if n is None and fraction is None:
        n = min(50_000, adata.n_obs)
    if fraction is not None:
        n = max(1, int(adata.n_obs * fraction))
    assert n is not None

    if n >= adata.n_obs:
        return adata.copy()

    groups = adata.obs[key]
    indices: list[int] = []
    for val in groups.unique():
        group_idx = np.where(groups == val)[0]
        group_n = max(1, int(len(group_idx) * n / adata.n_obs))
        group_n = min(group_n, len(group_idx))
        chosen = rng.choice(group_idx, size=group_n, replace=False)
        indices.extend(chosen.tolist())

    # Trim to exactly n cells randomly (not by index) to avoid truncation bias (BM-04)
    if len(indices) > n:
        indices = rng.choice(indices, size=n, replace=False).tolist()
    logger.info("Subsampled %d -> %d cells (stratified by %s)", adata.n_obs, len(indices), key)
    return adata[indices].copy()


def _resolve_best_method(comparison_df: pd.DataFrame) -> str:
    """Pick the best method from comparison results.

    Uses batch_score as primary metric (per Architecture.md Phase 3).
    Falls back to overall_score if batch_score not available.
    """
    if "batch_score" in comparison_df.columns:
        best_idx = comparison_df["batch_score"].idxmax()
    else:
        best_idx = comparison_df["overall_score"].idxmax()
    return str(comparison_df.loc[best_idx, "method"])


def _apply_integration_method(
    adata: AnnData,
    method: str,
    batch_key: str,
    celltype_key: str | None = None,
    use_gpu: str | bool = "auto",
    max_epochs: int = 200,
) -> str:
    """Apply a single integration method to adata. Returns obsm key."""
    # Resolve display name to short name
    short = _DISPLAY_TO_METHOD.get(method, method.lower())
    if short not in _METHOD_APPLY:
        raise ValueError(f"Unknown method: {method}. Valid: {list(_METHOD_APPLY)}")

    func_name, obsm_key = _METHOD_APPLY[short]
    if func_name is None:
        # PCA baseline -- nothing to do
        return obsm_key

    import importlib

    integrate_mod = importlib.import_module("sc_tools.pp.integrate")
    func = getattr(integrate_mod, func_name)

    kwargs: dict = {"batch_key": batch_key}
    if short in ("scvi", "scanvi", "cytovi", "resolvi"):
        kwargs["max_epochs"] = max_epochs
        kwargs["use_gpu"] = use_gpu
    if short == "scanvi" and celltype_key:
        kwargs["labels_key"] = celltype_key

    func(adata, **kwargs)
    return obsm_key


def run_full_integration_workflow(
    adata: AnnData,
    modality: str = "visium",
    batch_key: str = "library_id",
    celltype_key: str | None = None,
    methods: list[str] | None = None,
    output_dir: str | Path = "results",
    *,
    subsample_n: int | None = None,
    subsample_fraction: float | None = None,
    use_gpu: str | bool = "auto",
    max_epochs: int = 200,
    use_scib: str = "auto",
    save_intermediates: bool = True,
) -> tuple[AnnData, pd.DataFrame, str]:
    """Run the full integration benchmark workflow.

    Orchestrates: subsample -> benchmark all methods -> save intermediates
    -> select best (by batch score) -> apply to full dataset.

    Parameters
    ----------
    adata
        AnnData with raw or normalized data. Modified in place with the
        winning integration embedding.
    modality
        Data modality (determines default method set).
    batch_key
        Batch column in ``obs``.
    celltype_key
        Cell type column (optional; improves scoring but not required).
    methods
        Integration methods to benchmark. If ``None``, uses modality defaults.
    output_dir
        Directory for intermediate outputs.
    subsample_n
        Number of cells to subsample for benchmarking. Default: auto
        (all if <50k cells, else 50k stratified by batch).
    subsample_fraction
        Fraction of cells to subsample (overrides subsample_n).
    use_gpu
        GPU setting for VAE methods.
    max_epochs
        Max training epochs for VAE methods.
    use_scib
        Metric computation backend.
    save_intermediates
        If True, save per-method embeddings to
        ``output_dir/tmp/integration_test/{method}.h5ad``.

    Returns
    -------
    tuple[AnnData, pd.DataFrame, str]
        The AnnData with the best integration applied, the comparison
        DataFrame, and the name of the selected method.
    """
    output_dir = Path(output_dir)

    # Step 1: Subsample for benchmark
    if subsample_n is not None or subsample_fraction is not None or adata.n_obs > 50_000:
        subsample = _stratified_subsample(
            adata, batch_key, n=subsample_n, fraction=subsample_fraction
        )
    else:
        subsample = adata.copy()

    # Step 2: Run benchmark on subsample
    subsample, comparison_df = run_integration_benchmark(
        subsample,
        modality=modality,
        batch_key=batch_key,
        celltype_key=celltype_key,
        methods=methods,
        use_gpu=use_gpu,
        max_epochs=max_epochs,
        use_scib=use_scib,
    )

    if comparison_df.empty:
        raise RuntimeError("Integration benchmark produced no results")

    # Step 3: Save intermediates
    if save_intermediates:
        test_dir = output_dir / "tmp" / "integration_test"
        test_dir.mkdir(parents=True, exist_ok=True)

        for _, row in comparison_df.iterrows():
            method_name = str(row["method"])
            emb_key = str(row.get("embedding_key", ""))
            if emb_key and emb_key in subsample.obsm:
                try:
                    method_adata = subsample.copy()
                    method_adata.write_h5ad(test_dir / f"{method_name}.h5ad")
                    logger.info("Saved intermediate: %s", test_dir / f"{method_name}.h5ad")
                except Exception:
                    logger.warning("Failed to save intermediate for %s", method_name, exc_info=True)

    # Step 4: Select best method (batch score primary)
    best_method = _resolve_best_method(comparison_df)
    logger.info("Selected best integration method: %s", best_method)

    # Step 5: Record selection
    method_file = output_dir / "integration_method.txt"
    method_file.parent.mkdir(parents=True, exist_ok=True)
    method_file.write_text(best_method)

    # Step 6: Apply best method to full dataset
    try:
        obsm_key = _apply_integration_method(
            adata,
            best_method,
            batch_key=batch_key,
            celltype_key=celltype_key,
            use_gpu=use_gpu,
            max_epochs=max_epochs,
        )
        logger.info("Applied %s to full dataset -> %s", best_method, obsm_key)
    except Exception:
        logger.warning(
            "Failed to apply %s to full dataset; subsample results available",
            best_method,
            exc_info=True,
        )

    return adata, comparison_df, best_method
