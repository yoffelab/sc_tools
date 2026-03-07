"""
Main dispatcher for automated cell typing.

Usage
-----
>>> import sc_tools as st
>>> st.tl.annotate_celltypes(adata, method="sctype", marker_db=my_db)
>>> st.tl.annotate_celltypes(adata, method="celltypist", model="Immune_All_Low.pkl")
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import anndata as ad

from ._base import _store_results, get_backend


def annotate_celltypes(
    adata: ad.AnnData,
    method: str,
    *,
    cluster_key: str = "leiden",
    result_key: str | None = None,
    store_proba: bool = False,
    copy: bool = False,
    **kwargs,
) -> ad.AnnData:
    """
    Dispatch automated cell-type annotation to the requested backend.

    Results are written to:
    - ``adata.obs[f'celltype_auto_{key}']`` -- Categorical labels
    - ``adata.obs[f'celltype_auto_{key}_score']`` -- float64 confidence
    - ``adata.obsm[f'celltype_proba_{key}']`` -- probabilities (if store_proba)
    - ``adata.uns[f'celltype_auto_{key}']`` -- method metadata

    Parameters
    ----------
    adata
        Annotated data matrix.
    method
        Backend name (e.g. ``'sctype'``, ``'celltypist'``, ``'custom_gates'``, ``'ensemble'``).
    cluster_key
        Column in ``adata.obs`` with cluster labels (default: ``'leiden'``).
    result_key
        Override the suffix used for result keys. Defaults to ``method``.
    store_proba
        Whether to store per-cell probability DataFrame in ``obsm``.
    copy
        Return a copy of adata rather than modifying in-place.
    **kwargs
        Backend-specific keyword arguments passed through unchanged.

    Returns
    -------
    AnnData
        Annotated adata (same object if copy=False, new object if copy=True).
    """
    if copy:
        adata = adata.copy()

    backend_cls = get_backend(method)
    labels, scores, proba, metadata = backend_cls.run(
        adata,
        cluster_key=cluster_key,
        store_proba=store_proba,
        **kwargs,
    )

    key = result_key if result_key is not None else method
    _store_results(adata, key, labels, scores, proba, metadata, store_proba=store_proba)

    return adata
