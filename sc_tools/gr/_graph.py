"""
sc_tools.gr._graph — spatial_neighbors wrapper.

Thin wrapper around sq.gr.spatial_neighbors that:
- Always passes library_key
- Infers coord_type from adata.uns key presence
- Applies modality-specific defaults
- Stores params in adata.uns['gr']['spatial_neighbors']['params']
"""

from __future__ import annotations

import anndata as ad

try:
    import squidpy as sq

    SQUIDPY_AVAILABLE = True
except ImportError:
    SQUIDPY_AVAILABLE = False
    sq = None


def spatial_neighbors(
    adata: ad.AnnData,
    library_key: str = "library_id",
    coord_type: str | None = None,
    n_rings: int | None = None,
    delaunay: bool | None = None,
    **kwargs,
) -> None:
    """
    Compute spatial neighbors with per-ROI awareness.

    Thin wrapper around sq.gr.spatial_neighbors. By default:
    - Visium (coord_type grid detected from adata.uns): n_rings=1
    - All others: delaunay=True, coord_type=generic

    Results are stored in adata.obsp['spatial_connectivities'] (block-diagonal
    across ROIs when library_key is passed).

    Parameters
    ----------
    adata
        Multi-ROI AnnData with obsm['spatial'] and obs[library_key].
    library_key
        obs column identifying each ROI. Passed to sq.gr.spatial_neighbors
        to prevent cross-ROI edges.
    coord_type
        Override automatic coordinate type detection.
    n_rings
        Number of grid rings (Visium default: 1).
    delaunay
        Use Delaunay triangulation (non-Visium default: True).
    **kwargs
        Additional keyword arguments passed to sq.gr.spatial_neighbors.
    """
    if not SQUIDPY_AVAILABLE:
        raise ImportError("squidpy is required for sc_tools.gr.spatial_neighbors")

    # Modality-specific defaults
    modality = adata.uns.get("modality", "")
    is_visium = modality in ("visium", "visium_hd", "visium_hd_cell") or (
        "spatial" in adata.uns
        and any(
            "tissue_hires_scalef" in (v.get("scalefactors", {}) if isinstance(v, dict) else {})
            for v in adata.uns.get("spatial", {}).values()
            if isinstance(v, dict)
        )
    )

    if coord_type is None:
        coord_type = "grid" if is_visium else "generic"

    sq_kwargs: dict = {
        "library_key": library_key,
        "coord_type": coord_type,
    }

    if coord_type == "grid":
        sq_kwargs["n_rings"] = n_rings if n_rings is not None else 1
    else:
        sq_kwargs["delaunay"] = delaunay if delaunay is not None else True

    sq_kwargs.update(kwargs)

    # squidpy requires library_key to be categorical
    if library_key in adata.obs.columns and not hasattr(adata.obs[library_key], "cat"):
        adata.obs[library_key] = adata.obs[library_key].astype("category")

    sq.gr.spatial_neighbors(adata, **sq_kwargs)

    # Store params (build cleanly without redundancy)
    stored_params: dict = {
        "library_key": library_key,
        "coord_type": coord_type,
    }
    if coord_type == "grid":
        stored_params["n_rings"] = sq_kwargs.get("n_rings", 1)
    else:
        stored_params["delaunay"] = sq_kwargs.get("delaunay", True)
    # Include any caller-supplied extra kwargs (but not the already-stored keys)
    for k, v in kwargs.items():
        if k not in stored_params:
            stored_params[k] = v
    adata.uns.setdefault("gr", {}).setdefault("spatial_neighbors", {})["params"] = stored_params
