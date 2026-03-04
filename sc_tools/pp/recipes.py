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

from .integrate import run_cytovi, run_harmony, run_scvi
from .normalize import (
    arcsinh_transform,
    backup_raw,
    filter_genes_by_pattern,
    log_transform,
    normalize_total,
    scale,
)
from .reduce import cluster, pca, run_utag

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
            **kwargs,
        )
    elif modality in ("xenium", "visium_hd_cell"):
        _recipe_xenium(
            adata,
            batch_key=batch_key,
            integration=integration,
            n_top_genes=n_top_genes,
            resolution=resolution,
            filter_patterns=filter_patterns,
            use_gpu=use_gpu,
            **kwargs,
        )
    elif modality == "cosmx":
        _recipe_cosmx(
            adata,
            batch_key=batch_key,
            integration=integration,
            n_top_genes=n_top_genes,
            resolution=resolution,
            filter_patterns=filter_patterns,
            use_gpu=use_gpu,
            **kwargs,
        )
    elif modality == "imc":
        _recipe_imc(
            adata,
            batch_key=batch_key,
            integration=integration,
            resolution=resolution,
            use_gpu=use_gpu,
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
    **kwargs: Any,
) -> None:
    """Visium / Visium HD preprocessing recipe.

    1. Backup raw
    2. Filter MT/RP/HB genes
    3. HVG selection (with batch_key)
    4. Subset to HVGs
    5. scVI integration (from raw counts)
    6. PCA (for non-scVI visualizations)
    7. Neighbors + Leiden + UMAP
    """
    import scanpy as sc

    backup_raw(adata)
    filter_genes_by_pattern(adata, patterns=filter_patterns)

    hvg_batch_key = kwargs.get("hvg_batch_key", batch_key)

    if integration == "scvi":
        # scVI uses raw counts; HVG with seurat_v3 works on raw counts
        hvg_flavor = kwargs.get("hvg_flavor", "seurat_v3")
        sc.pp.highly_variable_genes(
            adata,
            flavor=hvg_flavor,
            n_top_genes=n_top_genes,
            batch_key=hvg_batch_key,
        )
        adata._inplace_subset_var(adata.var["highly_variable"].values)
        logger.info("Subset to %d HVGs", adata.n_vars)

        scvi_kwargs = {
            k: kwargs[k]
            for k in ("n_latent", "n_layers", "n_hidden", "max_epochs", "early_stopping")
            if k in kwargs
        }
        run_scvi(adata, batch_key=batch_key, use_gpu=use_gpu, **scvi_kwargs)

        # PCA on normalized data for visualization (supplementary to scVI)
        adata_norm = adata.copy()
        sc.pp.normalize_total(adata_norm)
        sc.pp.log1p(adata_norm)
        sc.tl.pca(adata_norm)
        adata.obsm["X_pca"] = adata_norm.obsm["X_pca"]
        adata.varm["PCs"] = adata_norm.varm["PCs"]
        adata.uns["pca"] = adata_norm.uns["pca"]
        del adata_norm
    else:
        # Non-scVI: normalize first, then HVG on log-normalized data
        normalize_total(adata)
        log_transform(adata)

        hvg_flavor = kwargs.get("hvg_flavor", "seurat")
        sc.pp.highly_variable_genes(
            adata,
            flavor=hvg_flavor,
            n_top_genes=n_top_genes,
            batch_key=hvg_batch_key,
        )
        adata._inplace_subset_var(adata.var["highly_variable"].values)
        logger.info("Subset to %d HVGs", adata.n_vars)

        n_comps = kwargs.get("n_comps", 50)
        pca(adata, n_comps=n_comps)

        if integration == "harmony":
            run_harmony(adata, batch_key=batch_key)

    # Clustering
    n_neighbors = kwargs.get("n_neighbors", 20)
    cluster(adata, resolution=resolution, n_neighbors=n_neighbors)


def _recipe_xenium(
    adata: AnnData,
    batch_key: str,
    integration: str,
    n_top_genes: int,
    resolution: float,
    filter_patterns: list[str] | None,
    use_gpu: str | bool,
    **kwargs: Any,
) -> None:
    """Xenium preprocessing recipe.

    1. Backup raw
    2. Filter unwanted genes
    3. Library-size normalization
    4. Skip log1p by default (per 2025 benchmarks; overridable)
    5. HVG selection
    6. Subset to HVGs
    7. Scale
    8. PCA
    9. Optional integration (Harmony or scVI)
    10. Neighbors + Leiden + UMAP
    """
    import scanpy as sc

    backup_raw(adata)
    filter_genes_by_pattern(adata, patterns=filter_patterns)

    normalize_total(adata)
    log_transform(adata)

    # HVG on log-normalized data
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    adata._inplace_subset_var(adata.var["highly_variable"].values)
    logger.info("Subset to %d HVGs", adata.n_vars)

    scale(adata, max_value=10)

    n_comps = kwargs.get("n_comps", 50)
    pca(adata, n_comps=n_comps)

    if integration == "harmony":
        run_harmony(adata, batch_key=batch_key)
    elif integration == "scvi":
        logger.warning(
            "scVI integration with Xenium: scVI expects raw counts. "
            "Consider using integration='harmony' or 'none' for Xenium."
        )
        scvi_kwargs = {
            k: kwargs[k]
            for k in ("n_latent", "n_layers", "n_hidden", "max_epochs", "early_stopping")
            if k in kwargs
        }
        run_scvi(adata, batch_key=batch_key, use_gpu=use_gpu, **scvi_kwargs)

    n_neighbors = kwargs.get("n_neighbors", 20)
    cluster(adata, resolution=resolution, n_neighbors=n_neighbors)


def _recipe_cosmx(
    adata: AnnData,
    batch_key: str,
    integration: str,
    n_top_genes: int,
    resolution: float,
    filter_patterns: list[str] | None,
    use_gpu: str | bool,
    **kwargs: Any,
) -> None:
    """CosMx preprocessing recipe.

    1. Backup raw
    2. Filter unwanted genes
    3. Library-size normalization
    4. Log1p
    5. HVG selection
    6. Subset to HVGs
    7. Scale
    8. PCA
    9. Optional integration (Harmony or scVI)
    10. Neighbors + Leiden + UMAP
    """
    import scanpy as sc

    backup_raw(adata)
    filter_genes_by_pattern(adata, patterns=filter_patterns)

    normalize_total(adata)
    log_transform(adata)

    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    adata._inplace_subset_var(adata.var["highly_variable"].values)
    logger.info("Subset to %d HVGs", adata.n_vars)

    scale(adata, max_value=10)

    n_comps = kwargs.get("n_comps", 50)
    pca(adata, n_comps=n_comps)

    if integration == "harmony":
        run_harmony(adata, batch_key=batch_key)
    elif integration == "scvi":
        logger.warning(
            "scVI integration with CosMx: scVI expects raw counts. "
            "Consider using integration='harmony' or 'none' for CosMx."
        )
        scvi_kwargs = {
            k: kwargs[k]
            for k in ("n_latent", "n_layers", "n_hidden", "max_epochs", "early_stopping")
            if k in kwargs
        }
        run_scvi(adata, batch_key=batch_key, use_gpu=use_gpu, **scvi_kwargs)

    n_neighbors = kwargs.get("n_neighbors", 20)
    cluster(adata, resolution=resolution, n_neighbors=n_neighbors)


def _recipe_imc(
    adata: AnnData,
    batch_key: str,
    integration: str,
    resolution: float,
    use_gpu: str | bool,
    **kwargs: Any,
) -> None:
    """IMC (Imaging Mass Cytometry) preprocessing recipe.

    1. Backup raw
    2. Arcsinh transform (cofactor=5) -- NOT log1p
    3. Scale
    4. PCA (n_comps capped to n_markers - 1)
    5. Optional integration (CytoVI or scVI fallback)
    6. Neighbors + Leiden + UMAP
    """
    backup_raw(adata)

    cofactor = kwargs.get("cofactor", 5)
    arcsinh_transform(adata, cofactor=cofactor)

    scale(adata, max_value=10)

    n_comps = kwargs.get("n_comps", min(20, adata.n_vars - 1))
    pca(adata, n_comps=n_comps, use_highly_variable=False)

    if integration == "cytovi":
        cytovi_kwargs = {k: kwargs[k] for k in ("n_latent", "max_epochs") if k in kwargs}
        run_cytovi(adata, batch_key=batch_key, use_gpu=use_gpu, **cytovi_kwargs)
    elif integration == "scvi":
        scvi_kwargs = {
            k: kwargs[k]
            for k in ("n_latent", "n_layers", "n_hidden", "max_epochs", "early_stopping")
            if k in kwargs
        }
        run_scvi(adata, batch_key=batch_key, use_gpu=use_gpu, **scvi_kwargs)
    elif integration == "harmony":
        run_harmony(adata, batch_key=batch_key)

    n_neighbors = kwargs.get("n_neighbors", 20)
    cluster(adata, resolution=resolution, n_neighbors=n_neighbors)
