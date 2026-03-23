"""Dimensionality reduction and clustering for preprocessing.

Provides:
- pca: Principal component analysis.
- neighbors: K-nearest neighbors graph.
- umap: UMAP embedding.
- leiden: Leiden clustering.
- cluster: Convenience wrapper (neighbors + leiden + umap).
- run_utag: Spatial-aware clustering via UTAG (soft dependency).
"""

from __future__ import annotations

import logging
from typing import Any

from anndata import AnnData

from ._gpu import get_backend

logger = logging.getLogger(__name__)

__all__ = [
    "pca",
    "neighbors",
    "umap",
    "leiden",
    "cluster",
    "run_utag",
]


def _auto_use_rep(adata: AnnData, use_rep: str | None) -> str:
    """Auto-detect the best representation for neighbor graph construction."""
    if use_rep is not None:
        return use_rep
    for key in ("X_scVI", "X_cytovi", "X_pca_harmony", "X_pca"):
        if key in adata.obsm:
            logger.info("Auto-detected use_rep='%s'", key)
            return key
    return "X_pca"


def pca(
    adata: AnnData,
    n_comps: int = 50,
    use_highly_variable: bool = True,
    **kwargs: Any,
) -> None:
    """Run PCA on the data.

    Wraps ``scanpy.tl.pca`` (or ``rapids_singlecell.tl.pca`` on GPU).

    Parameters
    ----------
    adata
        Annotated data. Modified in place.
    n_comps
        Number of principal components.
    use_highly_variable
        If True and ``adata.var['highly_variable']`` exists, use only HVGs.
    **kwargs
        Passed to the backend ``pca``.
    """
    backend, name = get_backend()

    # Clamp n_comps to valid range
    max_comps = min(adata.n_obs, adata.n_vars) - 1
    if n_comps > max_comps:
        logger.warning("Clamping n_comps from %d to %d (data shape)", n_comps, max_comps)
        n_comps = max(1, max_comps)

    hvg_available = "highly_variable" in adata.var.columns
    use_hvg = use_highly_variable and hvg_available

    logger.info("PCA (n_comps=%d, use_hvg=%s, backend=%s)", n_comps, use_hvg, name)
    backend.tl.pca(adata, n_comps=n_comps, use_highly_variable=use_hvg, **kwargs)


def neighbors(
    adata: AnnData,
    n_neighbors: int = 20,
    use_rep: str | None = None,
    **kwargs: Any,
) -> None:
    """Compute K-nearest neighbors graph.

    Wraps ``scanpy.pp.neighbors`` (or ``rapids_singlecell.pp.neighbors`` on GPU).
    Auto-detects ``use_rep``: X_scVI > X_cytovi > X_pca_harmony > X_pca.

    Parameters
    ----------
    adata
        Annotated data. Modified in place.
    n_neighbors
        Number of neighbors.
    use_rep
        Representation to use. Auto-detected if None.
    **kwargs
        Passed to the backend ``neighbors``.
    """
    use_rep = _auto_use_rep(adata, use_rep)
    backend, name = get_backend()
    logger.info("neighbors (n_neighbors=%d, use_rep='%s', backend=%s)", n_neighbors, use_rep, name)
    backend.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep, **kwargs)


def umap(adata: AnnData, **kwargs: Any) -> None:
    """Compute UMAP embedding.

    Wraps ``scanpy.tl.umap`` (or ``rapids_singlecell.tl.umap`` on GPU).

    Parameters
    ----------
    adata
        Annotated data with neighbor graph computed. Modified in place.
    **kwargs
        Passed to the backend ``umap``.
    """
    backend, name = get_backend()
    logger.info("UMAP (backend=%s)", name)
    backend.tl.umap(adata, **kwargs)


def leiden(
    adata: AnnData,
    resolution: float = 0.8,
    key_added: str = "leiden",
    **kwargs: Any,
) -> None:
    """Run Leiden clustering.

    Wraps ``scanpy.tl.leiden`` (or ``rapids_singlecell.tl.leiden`` on GPU).

    Parameters
    ----------
    adata
        Annotated data with neighbor graph computed. Modified in place.
    resolution
        Clustering resolution. Higher values yield more clusters.
    key_added
        Key in ``adata.obs`` to store cluster labels.
    **kwargs
        Passed to the backend ``leiden``.
    """
    backend, name = get_backend()
    logger.info(
        "Leiden clustering (resolution=%s, key='%s', backend=%s)", resolution, key_added, name
    )
    backend.tl.leiden(adata, resolution=resolution, key_added=key_added, **kwargs)


def cluster(
    adata: AnnData,
    resolution: float = 0.8,
    use_rep: str | None = None,
    n_neighbors: int = 20,
    key_added: str = "leiden",
    run_umap: bool = True,
    random_state: int = 0,
    **kwargs: Any,
) -> None:
    """Convenience: neighbors + leiden + umap in one call.

    Parameters
    ----------
    adata
        Annotated data. Modified in place.
    resolution
        Leiden resolution.
    use_rep
        Representation for neighbor graph. Auto-detected if None.
    n_neighbors
        Number of neighbors.
    key_added
        Key for cluster labels in ``adata.obs``.
    run_umap
        If True (default), compute UMAP embedding.
    random_state
        Random state for Leiden reproducibility (D-14, PRV-05).
    **kwargs
        Extra kwargs passed to ``neighbors()``.
    """
    neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep, **kwargs)
    leiden(adata, resolution=resolution, key_added=key_added, random_state=random_state)
    if run_umap:
        umap(adata)


def run_utag(
    adata: AnnData,
    max_dist: float = 20,
    slide_key: str | None = "library_id",
    clustering_method: str = "leiden",
    resolutions: list[float] | None = None,
    key_added: str = "utag",
    **kwargs: Any,
) -> None:
    """Run UTAG spatial-aware clustering.

    UTAG builds an adjacency graph from spatial coordinates, performs message
    passing to aggregate neighborhood features, then clusters to identify
    microanatomical domains. Runs *after* standard Leiden as a complementary
    spatial-aware annotation.

    Requires ``utag`` package:
    ``pip install git+https://github.com/ElementoLab/utag.git@main``

    Parameters
    ----------
    adata
        Annotated data with ``obsm['spatial']``. Modified in place.
    max_dist
        Distance threshold for cell adjacency. 10-20 for IMC, 10-100 for transcriptomics.
    slide_key
        Batch key for multi-image processing. None for single image.
    clustering_method
        Clustering method ("leiden" or "parc").
    resolutions
        List of resolutions to explore. Defaults to ``[0.5, 0.8, 1.0]``.
    key_added
        Prefix for cluster labels in ``adata.obs``.
    **kwargs
        Passed to ``utag.utag``.
    """
    try:
        import utag as utag_pkg
    except ImportError:
        raise ImportError(
            "UTAG is required for spatial clustering. Install with:\n"
            "  pip install git+https://github.com/ElementoLab/utag.git@main"
        ) from None

    if resolutions is None:
        resolutions = [0.5, 0.8, 1.0]

    logger.info(
        "UTAG (max_dist=%s, slide_key=%s, resolutions=%s)",
        max_dist,
        slide_key,
        resolutions,
    )

    utag_result = utag_pkg.utag(
        adata,
        max_dist=max_dist,
        slide_key=slide_key,
        clustering_method=clustering_method,
        resolutions=resolutions,
        **kwargs,
    )

    # Transfer UTAG cluster labels back to adata
    for col in utag_result.obs.columns:
        if col.startswith(clustering_method):
            new_key = f"{key_added}_{col}"
            adata.obs[new_key] = utag_result.obs[col].values
            logger.info("Added UTAG cluster labels: obs['%s']", new_key)
