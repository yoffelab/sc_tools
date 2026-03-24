"""Scale strategy pattern for preprocessing recipes."""

from __future__ import annotations

import logging
import warnings
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any

import numpy as np
from anndata import AnnData
from scipy.sparse import csr_matrix, issparse

from .integrate import run_cytovi, run_harmony, run_scvi
from .normalize import backup_raw, filter_genes_by_pattern, scale
from .projection import (
    _store_projection_info,
    build_knn_index,
    project_labels,
    project_representation,
    project_umap,
    subsample_stratified,
)
from .reduce import _auto_use_rep, leiden, neighbors, pca, umap

logger = logging.getLogger(__name__)


@dataclass
class SubsampleContext:
    """Mutable state flowing between strategy phases."""

    subsample_idx: np.ndarray | None = None
    knn_index: object | None = None
    use_rep: str = "X_pca"
    subsample_n: int = 1_000_000
    projection_k: int = 30
    random_state: int = 42
    stratify_key: str | None = None


class ScaleStrategy(ABC):
    """Defines how a preprocessing recipe executes based on data scale."""

    name: str

    @abstractmethod
    def prepare(
        self,
        adata: AnnData,
        *,
        raw_backup: str = "layer",
        filter_patterns: list[str] | None = None,
    ) -> None: ...

    @abstractmethod
    def select_features(
        self,
        adata: AnnData,
        *,
        n_top_genes: int = 2000,
        batch_key: str | None = None,
        **kw: Any,
    ) -> None: ...

    @abstractmethod
    def reduce_and_integrate(
        self,
        adata: AnnData,
        *,
        integration: str = "harmony",
        batch_key: str | None = None,
        n_comps: int = 50,
        **kw: Any,
    ) -> SubsampleContext | None: ...

    @abstractmethod
    def embed_and_cluster(
        self,
        adata: AnnData,
        *,
        ctx: SubsampleContext | None = None,
        resolution: float = 1.0,
        n_neighbors: int = 15,
        **kw: Any,
    ) -> None: ...


class SmallStrategy(ScaleStrategy):
    """Current code path -- all operations on full data in memory."""

    name = "small"

    def prepare(self, adata, *, raw_backup="layer", filter_patterns=None):
        backup_raw(adata)
        if filter_patterns is not None:
            filter_genes_by_pattern(adata, patterns=filter_patterns)

    def select_features(self, adata, *, n_top_genes=2000, batch_key=None, **kw):
        import scanpy as sc

        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key=batch_key, **kw)
        adata._inplace_subset_var(adata.var["highly_variable"].values)
        logger.info("Subset to %d HVGs", adata.n_vars)

    def reduce_and_integrate(
        self, adata, *, integration="harmony", batch_key=None, n_comps=50, **kw
    ):
        if integration == "scvi":
            scvi_kwargs = {
                k: kw[k]
                for k in (
                    "n_latent",
                    "n_layers",
                    "n_hidden",
                    "max_epochs",
                    "early_stopping",
                )
                if k in kw
            }
            use_gpu = kw.get("use_gpu", "auto")
            run_scvi(adata, batch_key=batch_key, use_gpu=use_gpu, **scvi_kwargs)
            # PCA on normalized copy for visualization
            import scanpy as sc

            adata_norm = adata.copy()
            sc.pp.normalize_total(adata_norm)
            sc.pp.log1p(adata_norm)
            sc.tl.pca(adata_norm, n_comps=n_comps)
            adata.obsm["X_pca"] = adata_norm.obsm["X_pca"]
            adata.varm["PCs"] = adata_norm.varm["PCs"]
            adata.uns["pca"] = adata_norm.uns["pca"]
            del adata_norm
        elif integration == "cytovi":
            cytovi_kwargs = {k: kw[k] for k in ("n_latent", "max_epochs") if k in kw}
            use_gpu = kw.get("use_gpu", "auto")
            run_cytovi(adata, batch_key=batch_key, use_gpu=use_gpu, **cytovi_kwargs)
        else:
            # none or harmony: scale + PCA
            scale(adata, max_value=10)
            # NaN cleanup for zero-variance features (common in IMC)
            _X = adata.X.toarray() if issparse(adata.X) else adata.X
            _nan_mask = np.any(np.isnan(_X), axis=0)
            if _nan_mask.any():
                warnings.warn(
                    f"reduce_and_integrate: {_nan_mask.sum()} zero-variance features "
                    "after scale; replacing NaN with 0.",
                    UserWarning,
                    stacklevel=2,
                )
                _X[:, _nan_mask] = 0.0
                adata.X = _X if not issparse(adata.X) else csr_matrix(_X)
            pca(adata, n_comps=n_comps)
            if integration == "harmony":
                run_harmony(adata, batch_key=batch_key)
        return None

    def embed_and_cluster(self, adata, *, ctx=None, resolution=1.0, n_neighbors=15, **kw):
        use_rep = _auto_use_rep(adata, use_rep=None)
        neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep)
        leiden(adata, resolution=resolution)
        umap(adata)


def _sparse_pca(adata, n_comps=50, gpu=False):
    """PCA on sparse data without densification."""
    if gpu:
        try:
            import cuml.decomposition
            from cupyx.scipy.sparse import csr_matrix as cu_csr

            X_gpu = cu_csr(adata.X)
            svd = cuml.decomposition.TruncatedSVD(n_components=n_comps)
            adata.obsm["X_pca"] = svd.fit_transform(X_gpu).get()
            logger.info("Sparse PCA via cuML TruncatedSVD (n_comps=%d)", n_comps)
            return
        except (ImportError, TypeError):
            logger.info("cuML sparse PCA unavailable, falling back to scanpy arpack")

    import scanpy as sc

    # Clamp n_comps to valid range (same guard as reduce.pca)
    max_comps = min(adata.n_obs, adata.n_vars) - 1
    if n_comps > max_comps:
        logger.warning("Clamping n_comps from %d to %d (data shape)", n_comps, max_comps)
        n_comps = max(1, max_comps)

    sc.pp.pca(adata, n_comps=n_comps, svd_solver="arpack", zero_center=True)
    logger.info("Sparse PCA via scanpy arpack solver (n_comps=%d)", n_comps)


class LargeStrategy(ScaleStrategy):
    """Subsample + project strategy for large datasets (7M+ cells)."""

    name = "large"

    def __init__(
        self,
        *,
        subsample_n=1_000_000,
        subsample_method="stratified",
        projection_k=30,
        random_state=42,
        has_gpu=False,
        backed=False,
    ):
        self.subsample_n = subsample_n
        self.subsample_method = subsample_method
        self.projection_k = projection_k
        self.random_state = random_state
        self.has_gpu = has_gpu
        self.backed = backed
        # Warn if FAISS unavailable
        try:
            import faiss  # noqa: F401
        except ImportError:
            warnings.warn(
                "FAISS not installed. kNN projection will fall back to sklearn "
                "BallTree, which is significantly slower for large datasets. "
                "Install with: pip install faiss-cpu (or faiss-gpu)",
                UserWarning,
                stacklevel=2,
            )

    def prepare(self, adata, *, raw_backup="layer", filter_patterns=None):
        # Layer-based backup: avoids doubling memory (no .raw)
        adata.layers["raw_counts"] = adata.X.copy()
        if filter_patterns is not None:
            filter_genes_by_pattern(adata, patterns=filter_patterns)

    def select_features(self, adata, *, n_top_genes=2000, batch_key=None, **kw):
        # Same as SmallStrategy -- HVG on full sparse is feasible
        import scanpy as sc

        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key=batch_key, **kw)
        adata._inplace_subset_var(adata.var["highly_variable"].values)
        logger.info("Subset to %d HVGs", adata.n_vars)

    def reduce_and_integrate(
        self, adata, *, integration="harmony", batch_key=None, n_comps=50, **kw
    ):
        if integration in ("scvi", "cytovi"):
            # scVI/CytoVI mini-batch naturally -- run on full data
            if integration == "cytovi":
                cytovi_kwargs = {k: kw[k] for k in ("n_latent", "max_epochs") if k in kw}
                use_gpu = kw.get("use_gpu", "auto")
                run_cytovi(adata, batch_key=batch_key, use_gpu=use_gpu, **cytovi_kwargs)
            else:
                scvi_kwargs = {
                    k: kw[k]
                    for k in (
                        "n_latent",
                        "n_layers",
                        "n_hidden",
                        "max_epochs",
                        "early_stopping",
                    )
                    if k in kw
                }
                use_gpu = kw.get("use_gpu", "auto")
                run_scvi(adata, batch_key=batch_key, use_gpu=use_gpu, **scvi_kwargs)
            # Skip supplementary PCA -- use latent space directly
            use_rep = _auto_use_rep(adata, use_rep=None)
            return SubsampleContext(
                subsample_n=self.subsample_n,
                projection_k=self.projection_k,
                random_state=self.random_state,
                stratify_key=batch_key,
                use_rep=use_rep,
            )

        # Sparse PCA (skip scale to avoid densification)
        _sparse_pca(adata, n_comps=n_comps, gpu=self.has_gpu)

        if integration == "none":
            return None  # No subsampling needed for no integration

        # Harmony: subsample -> correct -> project back
        ctx = SubsampleContext(
            subsample_n=self.subsample_n,
            projection_k=self.projection_k,
            random_state=self.random_state,
            stratify_key=batch_key,
        )
        ctx.subsample_idx = subsample_stratified(
            adata,
            n=ctx.subsample_n,
            stratify_key=batch_key,
            random_state=ctx.random_state,
        )
        adata_sub = adata[ctx.subsample_idx].copy()
        run_harmony(adata_sub, batch_key=batch_key)
        ctx.use_rep = _auto_use_rep(adata_sub, use_rep=None)
        ctx.knn_index = build_knn_index(adata_sub, use_rep="X_pca")
        project_representation(
            adata, ctx.knn_index, adata_sub, rep_key=ctx.use_rep, k=ctx.projection_k
        )
        return ctx

    def embed_and_cluster(self, adata, *, ctx=None, resolution=1.0, n_neighbors=15, **kw):
        # Build SubsampleContext if not provided
        if ctx is None:
            batch_key = kw.get("batch_key")
            ctx = SubsampleContext(
                subsample_n=self.subsample_n,
                projection_k=self.projection_k,
                random_state=self.random_state,
                stratify_key=batch_key,
                use_rep=_auto_use_rep(adata, use_rep=None),
            )

        # Subsample if not already done
        if ctx.subsample_idx is None:
            ctx.subsample_idx = subsample_stratified(
                adata,
                n=ctx.subsample_n,
                stratify_key=ctx.stratify_key,
                random_state=ctx.random_state,
            )

        adata_sub = adata[ctx.subsample_idx].copy()
        use_rep = ctx.use_rep if ctx.use_rep else _auto_use_rep(adata_sub, use_rep=None)
        neighbors(adata_sub, n_neighbors=n_neighbors, use_rep=use_rep)
        leiden(adata_sub, resolution=resolution)
        umap(adata_sub)

        # Build kNN index if not reused from reduce_and_integrate
        if ctx.knn_index is None:
            ctx.knn_index = build_knn_index(adata_sub, use_rep=use_rep)

        project_labels(adata, ctx.knn_index, adata_sub, label_key="leiden", k=ctx.projection_k)
        project_umap(adata, ctx.knn_index, adata_sub, k=ctx.projection_k)

        # Set leiden from projected labels
        adata.obs["leiden"] = adata.obs["leiden_projected"]

        _store_projection_info(adata, ctx, strategy_name=self.name)


# ---------------------------------------------------------------------------
# Auto-selection
# ---------------------------------------------------------------------------

# Platform-specific density estimates for memory calculation
PLATFORM_DENSITY = {
    # platform: (typical n_vars after HVG, sparsity, dtype_bytes)
    "visium_hd": (2000, 0.85, 4),
    "visium": (2000, 0.85, 4),
    "xenium": (300, 0.70, 4),
    "cosmx": (1000, 0.80, 4),
    "imc": (50, 0.0, 4),
    "merfish": (500, 0.75, 4),
}


def select_strategy(adata, config=None, platform=None):
    """Pick strategy based on estimated memory footprint + hardware + config.

    Parameters
    ----------
    adata
        AnnData object.
    config
        Optional dict with keys like ``backend`` ("auto"/"small"/"large"),
        ``large_strategy`` (dict of LargeStrategy kwargs).
    platform
        Platform name for density estimation.

    Returns
    -------
    ScaleStrategy instance.
    """
    from ._gpu import has_rapids

    # Config override takes priority
    if config and config.get("backend") not in (None, "auto"):
        return _build_strategy(config["backend"], config)

    # Use platform density table if available
    if platform and platform in PLATFORM_DENSITY:
        n_vars_est, _sparsity, dtype_bytes = PLATFORM_DENSITY[platform]
    else:
        n_vars_est = config.get("n_top_genes", 2000) if config else 2000
        dtype_bytes = 4

    # Peak memory: dense worst case for scale() in SmallStrategy
    dense_peak_gb = (adata.n_obs * n_vars_est * dtype_bytes) / 1e9
    available_gb = _get_available_memory_gb()
    has_gpu = has_rapids()

    # If dense matrix would exceed 50% of available memory -> LargeStrategy
    if dense_peak_gb > available_gb * 0.5:
        logger.info(
            "Estimated peak %.1fGB exceeds 50%% of %.0fGB available -> LargeStrategy",
            dense_peak_gb,
            available_gb,
        )
        return LargeStrategy(has_gpu=has_gpu, **_get_large_config(config))

    # Backed AnnData always routes to Large
    if adata.isbacked:
        return LargeStrategy(has_gpu=has_gpu, backed=True, **_get_large_config(config))

    return SmallStrategy()


def _get_available_memory_gb():
    """Return available system RAM in GB."""
    try:
        import psutil

        return psutil.virtual_memory().available / 1e9
    except ImportError:
        # Fallback: assume 16GB available
        logger.debug("psutil not available; assuming 16GB RAM")
        return 16.0


def _estimate_sparsity(X):
    """Return fraction of zeros in matrix."""
    if issparse(X):
        return 1.0 - (X.nnz / (X.shape[0] * X.shape[1]))
    return (X == 0).mean()


def _build_strategy(name, config=None):
    """Factory: map strategy name to instance."""
    if name == "small":
        return SmallStrategy()
    elif name == "large":
        from ._gpu import has_rapids

        return LargeStrategy(has_gpu=has_rapids(), **_get_large_config(config))
    else:
        raise ValueError(f"Unknown strategy backend '{name}'. Choose from: 'small', 'large'.")


def _get_large_config(config):
    """Extract LargeStrategy kwargs from config dict with defaults."""
    if not config:
        return {}
    large_cfg = config.get("large_strategy", {})
    result = {}
    for key in ("subsample_n", "subsample_method", "projection_k", "random_state"):
        if key in large_cfg:
            result[key] = large_cfg[key]
    return result
