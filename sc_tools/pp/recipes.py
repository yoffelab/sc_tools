"""Modality-aware preprocessing recipes.

Provides ``preprocess()`` as the main entry point, which dispatches to a
modality-specific recipe (visium, visium_hd, xenium, cosmx, imc). Each recipe
orchestrates the appropriate normalization, feature selection, integration,
dimensionality reduction, and clustering steps.

Individual steps are all importable from ``sc_tools.pp`` for fine-grained control.
"""

from __future__ import annotations

import logging
from typing import Any

from anndata import AnnData

from .normalize import (
    log_transform,
    normalize_total,
)
from .reduce import run_utag
from .strategy import SmallStrategy, select_strategy

logger = logging.getLogger(__name__)

__all__ = ["preprocess"]

VALID_MODALITIES = {"visium", "visium_hd", "visium_hd_cell", "xenium", "cosmx", "imc"}
VALID_INTEGRATIONS = {"scvi", "harmony", "cytovi", "none"}


def preprocess(
    adata: AnnData,
    modality: str = "visium",
    batch_key: str = "library_id",
    integration: str = "scvi",
    spatial_clustering: str | None = None,
    n_top_genes: int = 2000,
    resolution: float = 0.8,
    filter_patterns: list[str] | None = None,
    use_gpu: str | bool = "auto",
    copy: bool = False,
    **kwargs: Any,
) -> AnnData:
    """Modality-aware preprocessing pipeline.

    Dispatches to a recipe based on ``modality``. Each recipe handles
    normalization, feature selection, batch integration, dimensionality
    reduction, and clustering appropriate for the data type.

    Parameters
    ----------
    adata
        Annotated data with raw counts/intensities in ``X``.
    modality
        Data modality: ``"visium"``, ``"visium_hd"``, ``"xenium"``, ``"cosmx"``,
        or ``"imc"``.
    batch_key
        Column in ``adata.obs`` for batch/library identification.
    integration
        Integration method: ``"scvi"`` (default), ``"harmony"``, ``"cytovi"``
        (IMC only), or ``"none"``.
    spatial_clustering
        If ``"utag"``, run UTAG spatial-aware clustering after standard Leiden.
    n_top_genes
        Number of highly variable genes to select (not used for IMC).
    resolution
        Leiden clustering resolution.
    filter_patterns
        Gene name patterns to exclude. None uses modality defaults.
    use_gpu
        ``"auto"`` (detect), True, or False.
    copy
        If True, operate on a copy and return it.
    **kwargs
        Extra arguments passed to integration and recipe-specific steps.
        Recognized keys:

        - ``strategy``: explicit ``ScaleStrategy`` instance (default: auto-select).
        - ``n_latent``, ``n_layers``, ``n_hidden``, ``max_epochs``,
          ``early_stopping``: passed to ``run_scvi``/``run_cytovi``.
        - ``n_neighbors``: passed to ``neighbors()`` (default 20).
        - ``n_comps``: number of PCA components (default 50).
        - ``cofactor``: arcsinh cofactor for IMC (default 5).
        - ``skip_log1p``: if True, skip log1p in Xenium recipe (default True).
        - ``max_dist``, ``utag_resolutions``: passed to ``run_utag``.
        - ``hvg_flavor``, ``hvg_batch_key``: passed to HVG selection.

    Returns
    -------
    AnnData
        Preprocessed data with clustering, embeddings, and (optionally)
        batch-corrected latent space.
    """
    modality = modality.lower()
    integration = integration.lower()

    if modality not in VALID_MODALITIES:
        raise ValueError(f"Unknown modality '{modality}'. Choose from {VALID_MODALITIES}")
    if integration not in VALID_INTEGRATIONS:
        raise ValueError(f"Unknown integration '{integration}'. Choose from {VALID_INTEGRATIONS}")

    if copy:
        adata = adata.copy()

    # Extract strategy from kwargs (not a named param to keep backward compat)
    strategy = kwargs.pop("strategy", None)
    if strategy is None:
        strategy = select_strategy(adata, config=kwargs.pop("config", None), platform=modality)

    logger.info("Preprocessing: modality=%s, integration=%s", modality, integration)

    if modality in ("visium", "visium_hd"):
        _recipe_visium(
            adata,
            batch_key=batch_key,
            integration=integration,
            n_top_genes=n_top_genes,
            resolution=resolution,
            filter_patterns=filter_patterns,
            use_gpu=use_gpu,
            strategy=strategy,
            **kwargs,
        )
    elif modality in ("xenium", "visium_hd_cell", "cosmx"):
        _recipe_targeted_panel(
            adata,
            batch_key=batch_key,
            integration=integration,
            n_top_genes=n_top_genes,
            resolution=resolution,
            filter_patterns=filter_patterns,
            use_gpu=use_gpu,
            strategy=strategy,
            modality=modality,
            **kwargs,
        )
    elif modality == "imc":
        _recipe_imc(
            adata,
            batch_key=batch_key,
            integration=integration,
            resolution=resolution,
            use_gpu=use_gpu,
            strategy=strategy,
            **kwargs,
        )

    # Optional spatial clustering
    if spatial_clustering == "utag":
        max_dist = kwargs.get("max_dist", 15 if modality == "imc" else 20)
        utag_resolutions = kwargs.get("utag_resolutions", None)
        run_utag(
            adata,
            max_dist=max_dist,
            slide_key=batch_key,
            resolutions=utag_resolutions,
        )

    return adata


def _recipe_visium(
    adata: AnnData,
    batch_key: str,
    integration: str,
    n_top_genes: int,
    resolution: float,
    filter_patterns: list[str] | None,
    use_gpu: str | bool,
    strategy: SmallStrategy | None = None,
    **kwargs: Any,
) -> None:
    """Visium / Visium HD preprocessing recipe.

    1. Backup raw + filter genes (strategy.prepare)
    2. HVG selection (strategy.select_features)
    3. Integration / dimensionality reduction (strategy.reduce_and_integrate)
    4. Neighbors + Leiden + UMAP (strategy.embed_and_cluster)
    """
    if strategy is None:
        strategy = SmallStrategy()

    strategy.prepare(adata, filter_patterns=filter_patterns)

    hvg_batch_key = kwargs.get("hvg_batch_key", batch_key)

    if integration == "scvi":
        # scVI uses raw counts; HVG with seurat_v3
        hvg_flavor = kwargs.get("hvg_flavor", "seurat_v3")
        strategy.select_features(
            adata, n_top_genes=n_top_genes, flavor=hvg_flavor, batch_key=hvg_batch_key
        )
        ctx = strategy.reduce_and_integrate(
            adata,
            integration="scvi",
            batch_key=batch_key,
            n_comps=kwargs.get("n_comps", 50),
            use_gpu=use_gpu,
            **{
                k: kwargs[k]
                for k in ("n_latent", "n_layers", "n_hidden", "max_epochs", "early_stopping")
                if k in kwargs
            },
        )
    else:
        # Non-scVI: normalize first, then HVG
        normalize_total(adata)
        log_transform(adata)
        hvg_flavor = kwargs.get("hvg_flavor", "seurat")
        strategy.select_features(
            adata, n_top_genes=n_top_genes, flavor=hvg_flavor, batch_key=hvg_batch_key
        )
        ctx = strategy.reduce_and_integrate(
            adata,
            integration=integration,
            batch_key=batch_key,
            n_comps=kwargs.get("n_comps", 50),
        )

    n_neighbors = kwargs.get("n_neighbors", 20)
    strategy.embed_and_cluster(adata, ctx=ctx, resolution=resolution, n_neighbors=n_neighbors)


def _recipe_targeted_panel(
    adata: AnnData,
    batch_key: str,
    integration: str,
    n_top_genes: int,
    resolution: float,
    filter_patterns: list[str] | None,
    use_gpu: str | bool,
    strategy: SmallStrategy | None = None,
    modality: str = "xenium",
    **kwargs: Any,
) -> None:
    """Targeted panel recipe (Xenium, CosMx, Visium HD cell-seg)."""
    if strategy is None:
        strategy = SmallStrategy()

    strategy.prepare(adata, filter_patterns=filter_patterns)

    if integration == "scvi":
        # scVI needs raw counts -- skip normalize, use seurat_v3 for HVG selection
        strategy.select_features(adata, n_top_genes=n_top_genes, flavor="seurat_v3")
        scvi_kwargs = {
            k: kwargs[k]
            for k in ("n_latent", "n_layers", "n_hidden", "max_epochs", "early_stopping")
            if k in kwargs
        }
        ctx = strategy.reduce_and_integrate(
            adata,
            integration="scvi",
            batch_key=batch_key,
            n_comps=kwargs.get("n_comps", 50),
            use_gpu=use_gpu,
            **scvi_kwargs,
        )
    else:
        normalize_total(adata)
        log_transform(adata)
        strategy.select_features(adata, n_top_genes=n_top_genes)
        ctx = strategy.reduce_and_integrate(
            adata,
            integration=integration,
            batch_key=batch_key,
            n_comps=kwargs.get("n_comps", 50),
        )

    n_neighbors = kwargs.get("n_neighbors", 20)
    strategy.embed_and_cluster(adata, ctx=ctx, resolution=resolution, n_neighbors=n_neighbors)


def _recipe_imc(
    adata: AnnData,
    batch_key: str,
    integration: str,
    resolution: float,
    use_gpu: str | bool,
    strategy: SmallStrategy | None = None,
    **kwargs: Any,
) -> None:
    """IMC (Imaging Mass Cytometry) preprocessing recipe.

    1. Backup raw (strategy.prepare)
    2. Arcsinh transform (recipe owns normalization)
    3. Scale + PCA + integration (strategy.reduce_and_integrate)
    4. Neighbors + Leiden + UMAP (strategy.embed_and_cluster)
    """
    if strategy is None:
        strategy = SmallStrategy()

    strategy.prepare(adata, filter_patterns=None)

    # IMC normalization (recipe owns this, not strategy)
    from .normalize import arcsinh_transform

    cofactor = kwargs.get("cofactor", 5)
    arcsinh_transform(adata, cofactor=cofactor)

    # IMC: strategy handles scale + PCA + integration
    # NaN cleanup after scale is handled inside strategy.reduce_and_integrate
    n_comps = kwargs.get("n_comps", min(20, adata.n_vars - 1))

    if integration == "cytovi":
        ctx = strategy.reduce_and_integrate(
            adata,
            integration="cytovi",
            batch_key=batch_key,
            n_comps=n_comps,
            use_gpu=use_gpu,
            **{k: kwargs[k] for k in ("n_latent", "max_epochs") if k in kwargs},
        )
    elif integration == "scvi":
        scvi_kwargs = {
            k: kwargs[k]
            for k in ("n_latent", "n_layers", "n_hidden", "max_epochs", "early_stopping")
            if k in kwargs
        }
        ctx = strategy.reduce_and_integrate(
            adata,
            integration="scvi",
            batch_key=batch_key,
            n_comps=n_comps,
            use_gpu=use_gpu,
            **scvi_kwargs,
        )
    else:
        ctx = strategy.reduce_and_integrate(
            adata,
            integration=integration,
            batch_key=batch_key,
            n_comps=n_comps,
        )

    n_neighbors = kwargs.get("n_neighbors", 20)
    strategy.embed_and_cluster(adata, ctx=ctx, resolution=resolution, n_neighbors=n_neighbors)
