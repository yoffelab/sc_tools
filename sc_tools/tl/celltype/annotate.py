"""
Main dispatcher for automated cell typing.

Includes panel-aware dispatch: targeted panels (n_vars < 1000) are restricted
to validated methods (sctype, custom_gates) unless force_method=True.

Usage
-----
>>> import sc_tools as st
>>> st.tl.annotate_celltypes(adata, method="sctype", marker_db=my_db)
>>> st.tl.annotate_celltypes(adata, method="celltypist", model="Immune_All_Low.pkl")
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import anndata as ad

from sc_tools.errors import SCToolsDataError

from ._base import _store_results, get_backend

logger = logging.getLogger(__name__)

# Panel dispatch constants (SCI-04, D-11, D-12, D-13)
PANEL_VALIDATED_METHODS = {"sctype", "custom_gates"}
WHOLE_TRANSCRIPTOME_METHODS = {"celltypist", "scgpt", "geneformer", "scarches", "singler"}
PANEL_THRESHOLD = 1000


def annotate_celltypes(
    adata: ad.AnnData,
    method: str,
    *,
    cluster_key: str = "leiden",
    result_key: str | None = None,
    store_proba: bool = False,
    copy: bool = False,
    force_method: bool = False,
    **kwargs,
) -> ad.AnnData:
    """
    Dispatch automated cell-type annotation to the requested backend.

    Panel guard: when ``adata.n_vars < 1000`` (or ``adata.raw.n_vars`` if
    available), whole-transcriptome methods are blocked unless
    ``force_method=True``.

    Results are written to:
    - ``adata.obs[f'celltype_auto_{key}']`` -- Categorical labels
    - ``adata.obs[f'celltype_auto_{key}_score']`` -- float64 confidence
    - ``adata.obsm[f'celltype_proba_{key}']`` -- probabilities (if store_proba)
    - ``adata.uns[f'celltype_auto_{key}']`` -- method metadata
    - ``adata.uns['panel_dispatch']`` -- panel detection provenance

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
    force_method
        If True, allow whole-transcriptome methods on panel data with a warning.
    **kwargs
        Backend-specific keyword arguments passed through unchanged.

    Returns
    -------
    AnnData
        Annotated adata (same object if copy=False, new object if copy=True).

    Raises
    ------
    SCToolsDataError
        If a whole-transcriptome method is requested for panel data without
        ``force_method=True``.
    """
    if copy:
        adata = adata.copy()

    # Panel detection (Pitfall 3: use raw.n_vars if available to avoid HVG false positives)
    n_vars_check = adata.n_vars
    if adata.raw is not None:
        n_vars_check = adata.raw.n_vars

    panel_detected = n_vars_check < PANEL_THRESHOLD

    # Panel guard (D-11, D-12)
    if panel_detected and method in WHOLE_TRANSCRIPTOME_METHODS:
        if force_method:
            logger.warning(
                "Panel detected (n_vars=%d) but force_method used for '%s'. "
                "Results may be unreliable.",
                n_vars_check,
                method,
            )
        else:
            raise SCToolsDataError(
                f"Method '{method}' requires whole-transcriptome data but "
                f"panel detected (n_vars={n_vars_check}). "
                f"Use sctype or custom_gates, or pass force_method=True.",
                suggestion="Use --method sctype or --force-method",
            )

    # Store panel dispatch provenance (D-13)
    adata.uns["panel_dispatch"] = {
        "panel_detected": panel_detected,
        "n_vars": adata.n_vars,
        "restricted_methods": sorted(PANEL_VALIDATED_METHODS) if panel_detected else [],
    }

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
