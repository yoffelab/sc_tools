"""sc_tools.pp: Modality-aware preprocessing for spatial omics data.

Provides normalization, batch integration, dimensionality reduction,
clustering, and full preprocessing recipes for Visium, Visium HD,
Xenium, CosMx, and IMC data.

GPU acceleration via rapids-singlecell is auto-detected; falls back to scanpy.

Usage (recipe)::

    import sc_tools.pp as pp
    adata = pp.preprocess(adata, modality="visium", batch_key="library_id")

Usage (individual steps)::

    import sc_tools.pp as pp
    pp.backup_raw(adata)
    pp.filter_genes_by_pattern(adata)
    pp.normalize_total(adata)
    pp.log_transform(adata)
    pp.pca(adata)
    pp.cluster(adata, resolution=0.8)
"""

from .integrate import (
    run_bbknn,
    run_combat,
    run_cytovi,
    run_harmony,
    run_resolvi,
    run_scanorama,
    run_scanvi,
    run_scvi,
)
from .integration_configs import get_resolvi_ss_config, get_scanvi_config
from .normalize import (
    arcsinh_transform,
    backup_raw,
    filter_genes_by_pattern,
    log_transform,
    normalize_total,
    scale,
)
from .recipes import preprocess
from .reduce import cluster, leiden, neighbors, pca, run_utag, umap

__all__ = [
    # Recipes
    "preprocess",
    # Normalization
    "backup_raw",
    "normalize_total",
    "log_transform",
    "scale",
    "arcsinh_transform",
    "filter_genes_by_pattern",
    # Integration
    "run_scvi",
    "run_harmony",
    "run_cytovi",
    "run_combat",
    "run_bbknn",
    "run_scanorama",
    "run_scanvi",
    "run_resolvi",
    # Integration configs
    "get_scanvi_config",
    "get_resolvi_ss_config",
    # Dimensionality reduction & clustering
    "pca",
    "neighbors",
    "umap",
    "leiden",
    "cluster",
    "run_utag",
]
