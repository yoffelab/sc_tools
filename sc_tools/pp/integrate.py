"""Batch integration methods for preprocessing.

Provides:
- run_scvi: scVI variational inference (for sequencing-based spatial data).
- run_harmony: Harmony integration (fast, PCA-based).
- run_cytovi: CytoVI for mass cytometry / IMC protein data.

All methods are soft dependencies: ImportError with install hints if not available.
"""

from __future__ import annotations

import logging
from typing import Any

from anndata import AnnData

logger = logging.getLogger(__name__)

__all__ = [
    "run_scvi",
    "run_harmony",
    "run_cytovi",
]


def _detect_gpu(use_gpu: str | bool) -> bool:
    """Resolve 'auto' GPU setting to True/False."""
    if isinstance(use_gpu, bool):
        return use_gpu
    # auto-detect
    try:
        import torch

        return torch.cuda.is_available()
    except ImportError:
        return False


def run_scvi(
    adata: AnnData,
    batch_key: str = "library_id",
    n_latent: int = 30,
    n_layers: int = 2,
    n_hidden: int = 128,
    max_epochs: int = 400,
    early_stopping: bool = True,
    use_gpu: str | bool = "auto",
    layer: str | None = None,
    **kwargs: Any,
) -> None:
    """Run scVI batch integration.

    Trains an scVI model on raw counts and stores the latent representation
    in ``adata.obsm['X_scVI']``. Model parameters are recorded in
    ``adata.uns['scvi_params']``.

    Requires ``scvi-tools``: ``pip install sc-tools[deconvolution]``

    Parameters
    ----------
    adata
        Annotated data with raw counts in ``X`` (unnormalized).
        Modified in place.
    batch_key
        Column in ``adata.obs`` for batch correction.
    n_latent
        Dimensionality of the latent space.
    n_layers
        Number of hidden layers.
    n_hidden
        Number of units per hidden layer.
    max_epochs
        Maximum training epochs.
    early_stopping
        If True, stop training when validation loss plateaus.
    use_gpu
        ``"auto"`` (default), True, or False.
    layer
        Layer in ``adata.layers`` to use. None uses ``X``.
    **kwargs
        Passed to ``scvi.model.SCVI()``.
    """
    try:
        import scvi
    except ImportError:
        raise ImportError(
            "scvi-tools is required for scVI integration. Install with:\n"
            "  pip install sc-tools[deconvolution]"
        ) from None

    gpu = _detect_gpu(use_gpu)
    logger.info(
        "scVI (batch_key='%s', n_latent=%d, max_epochs=%d, gpu=%s)",
        batch_key,
        n_latent,
        max_epochs,
        gpu,
    )

    scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key, layer=layer)

    model = scvi.model.SCVI(
        adata,
        n_latent=n_latent,
        n_layers=n_layers,
        n_hidden=n_hidden,
        **kwargs,
    )

    train_kwargs: dict[str, Any] = {"max_epochs": max_epochs}
    if early_stopping:
        train_kwargs["early_stopping"] = True
    if gpu:
        train_kwargs["accelerator"] = "gpu"

    model.train(**train_kwargs)

    adata.obsm["X_scVI"] = model.get_latent_representation()
    adata.uns["scvi_params"] = {
        "batch_key": batch_key,
        "n_latent": n_latent,
        "n_layers": n_layers,
        "n_hidden": n_hidden,
        "max_epochs": max_epochs,
        "early_stopping": early_stopping,
    }
    logger.info("scVI latent stored in obsm['X_scVI'] (shape=%s)", adata.obsm["X_scVI"].shape)


def run_harmony(
    adata: AnnData,
    batch_key: str = "library_id",
    basis: str = "X_pca",
    key_added: str = "X_pca_harmony",
    **kwargs: Any,
) -> None:
    """Run Harmony integration on a PCA embedding.

    Wraps ``scanpy.external.pp.harmony_integrate``. Stores corrected embedding
    in ``adata.obsm[key_added]``.

    Requires ``harmonypy``: ``pip install harmonypy``

    Parameters
    ----------
    adata
        Annotated data with PCA computed (``obsm[basis]``). Modified in place.
    batch_key
        Column in ``adata.obs`` for batch correction.
    basis
        PCA embedding key in ``obsm``.
    key_added
        Key for the corrected embedding in ``obsm``.
    **kwargs
        Passed to ``harmony_integrate``.
    """
    try:
        import harmonypy  # noqa: F401
    except ImportError:
        raise ImportError(
            "harmonypy is required for Harmony integration. Install with:\n  pip install harmonypy"
        ) from None

    import scanpy as sc

    if basis not in adata.obsm:
        raise ValueError(f"'{basis}' not found in adata.obsm. Run PCA first.")

    logger.info("Harmony (batch_key='%s', basis='%s')", batch_key, basis)

    sc.external.pp.harmony_integrate(
        adata,
        key=batch_key,
        basis=basis,
        adjusted_basis=key_added,
        **kwargs,
    )
    logger.info("Harmony corrected embedding stored in obsm['%s']", key_added)


def run_cytovi(
    adata: AnnData,
    batch_key: str = "library_id",
    n_latent: int = 20,
    max_epochs: int = 300,
    use_gpu: str | bool = "auto",
    **kwargs: Any,
) -> None:
    """Run CytoVI for mass cytometry / IMC protein data integration.

    CytoVI is a totalVI-inspired model from scvi-tools designed for protein
    marker data. Stores latent representation in ``adata.obsm['X_cytovi']``.

    Requires ``scvi-tools``: ``pip install sc-tools[deconvolution]``

    Parameters
    ----------
    adata
        Annotated data with protein intensities. Modified in place.
    batch_key
        Column in ``adata.obs`` for batch correction.
    n_latent
        Dimensionality of the latent space.
    max_epochs
        Maximum training epochs.
    use_gpu
        ``"auto"`` (default), True, or False.
    **kwargs
        Passed to the model constructor.
    """
    try:
        import scvi
    except ImportError:
        raise ImportError(
            "scvi-tools is required for CytoVI integration. Install with:\n"
            "  pip install sc-tools[deconvolution]"
        ) from None

    gpu = _detect_gpu(use_gpu)
    logger.info(
        "CytoVI (batch_key='%s', n_latent=%d, max_epochs=%d, gpu=%s)",
        batch_key,
        n_latent,
        max_epochs,
        gpu,
    )

    # CytoVI uses totalVI-style setup for protein data
    # Fall back to SCVI if CytoVI is not available in the installed version
    model_cls = getattr(scvi.model, "CytoVI", None)
    if model_cls is None:
        logger.warning(
            "CytoVI model not found in scvi-tools %s; falling back to SCVI for batch integration",
            getattr(scvi, "__version__", "unknown"),
        )
        run_scvi(
            adata,
            batch_key=batch_key,
            n_latent=n_latent,
            max_epochs=max_epochs,
            use_gpu=use_gpu,
            **kwargs,
        )
        # Rename key if scVI was used as fallback
        if "X_scVI" in adata.obsm:
            adata.obsm["X_cytovi"] = adata.obsm.pop("X_scVI")
        return

    model_cls.setup_anndata(adata, batch_key=batch_key)
    model = model_cls(adata, n_latent=n_latent, **kwargs)

    train_kwargs: dict[str, Any] = {"max_epochs": max_epochs}
    if gpu:
        train_kwargs["accelerator"] = "gpu"

    model.train(**train_kwargs)

    adata.obsm["X_cytovi"] = model.get_latent_representation()
    adata.uns["cytovi_params"] = {
        "batch_key": batch_key,
        "n_latent": n_latent,
        "max_epochs": max_epochs,
    }
    logger.info("CytoVI latent stored in obsm['X_cytovi'] (shape=%s)", adata.obsm["X_cytovi"].shape)
