"""
Spatial QC: spatially variable genes (squidpy).

Mitochondrial and hemoglobin percentages per spot are computed in
sc_tools.qc.metrics.calculate_qc_metrics (pct_counts_mt, pct_counts_hb in adata.obs).
"""

from __future__ import annotations

from typing import Any

import pandas as pd

from anndata import AnnData

__all__ = ["spatially_variable_genes", "spatially_variable_genes_per_library"]


def spatially_variable_genes(
    adata: AnnData,
    *,
    mode: str = "moran",
    coord_key: str = "spatial",
    genes: list[str] | None = None,
    n_top_genes: int | None = None,
    threshold_i: float | None = None,
    copy_var: bool = True,
    n_perms: int = 100,
    n_jobs: int = 1,
    **kwargs: Any,
) -> pd.DataFrame | None:
    """
    Compute spatially variable genes using squidpy (Moran's I or Geary's C).

    Builds spatial neighbors if missing, runs spatial autocorrelation, and
    optionally writes results to adata.var (spatial_moran_i, spatial_pval,
    spatially_variable).

    Parameters
    ----------
    adata : AnnData
        Annotated data with adata.obsm[coord_key] (default 'spatial').
    mode : str
        'moran' or 'geary' (default 'moran').
    coord_key : str
        Key in adata.obsm for coordinates (default 'spatial').
    genes : list of str or None
        Genes to test; if None, uses adata.var_names (or highly_variable if present).
    n_top_genes : int or None
        If set, mark this many top genes by statistic as spatially_variable.
    threshold_i : float or None
        If set, mark genes with Moran's I >= threshold_i as spatially_variable.
    copy_var : bool
        If True, copy I and pval from uns to adata.var and add spatially_variable (default True).
    n_perms : int
        Permutations for p-value (default 100). Passed to squidpy.
    n_jobs : int
        Parallel jobs (default 1). Passed to squidpy.
    **kwargs
        Passed to squidpy.gr.spatial_autocorr.

    Returns
    -------
    DataFrame or None
        Autocorrelation results (index = genes) if copy_var is False or for inspection;
        otherwise None (results in adata.uns and adata.var).
    """
    try:
        import squidpy as sq
    except ImportError as e:
        raise ImportError("squidpy is required for spatially_variable_genes") from e

    if coord_key not in adata.obsm:
        raise ValueError(f"adata.obsm[{coord_key!r}] not found. Required for spatial neighbors.")

    # Build spatial graph if missing
    if "spatial_neighbors" not in adata.obsp:
        sq.gr.spatial_neighbors(adata, coord_system="generic", coord_key=coord_key)

    if genes is not None:
        genes = [g for g in genes if g in adata.var_names]
        if not genes:
            genes = adata.var_names.tolist()
    else:
        if "highly_variable" in adata.var.columns and adata.var["highly_variable"].any():
            genes = adata.var_names[adata.var["highly_variable"]].tolist()
        else:
            genes = adata.var_names.tolist()

    sq.gr.spatial_autocorr(
        adata,
        mode=mode,
        genes=genes,
        n_perms=n_perms,
        n_jobs=n_jobs,
        **kwargs,
    )

    # squidpy stores in adata.uns['moranI'] or 'gearyC' (DataFrame with I, pval_norm, etc.)
    key = "moranI" if mode == "moran" else "gearyC"
    if key not in adata.uns:
        return None

    result = adata.uns[key]
    if not isinstance(result, pd.DataFrame):
        return None

    stat_col = "I" if mode == "moran" else "C"
    pval_col = "pval_norm" if "pval_norm" in result.columns else result.columns[-1]

    if copy_var:
        adata.var["spatial_" + stat_col.lower()] = adata.var_names.map(
            result[stat_col].reindex(adata.var_names).to_dict()
        ).astype(float)
        adata.var["spatial_pval"] = adata.var_names.map(
            result[pval_col].reindex(adata.var_names).to_dict()
        ).astype(float)

        spatially_var = pd.Series(False, index=adata.var_names)
        if n_top_genes is not None and n_top_genes > 0:
            top = result.nlargest(n_top_genes, stat_col).index
            spatially_var.loc[spatially_var.index.isin(top)] = True
        if threshold_i is not None:
            col = "spatial_" + stat_col.lower()
            spatially_var = spatially_var | (adata.var[col] >= threshold_i)
        adata.var["spatially_variable"] = spatially_var.values

    return result


def spatially_variable_genes_per_library(
    adata: AnnData,
    library_id_col: str = "library_id",
    *,
    mode: str = "moran",
    coord_key: str = "spatial",
    n_top_genes: int | None = None,
    threshold_i: float | None = None,
    n_perms: int = 100,
    n_jobs: int = 1,
    **kwargs: Any,
) -> dict[str, pd.DataFrame] | None:
    """
    Run spatially_variable_genes on each library (sample) separately and store results in uns.

    Spatial neighbors are built per sample, so each library is subset and processed
    independently. If library_id_col is not in adata.obs, returns None without error
    (caller should skip SVG and continue other QC).

    Parameters
    ----------
    adata : AnnData
        Annotated data with adata.obs[library_id_col] and adata.obsm[coord_key].
    library_id_col : str
        Column in adata.obs identifying library/sample (default 'library_id').
    mode, coord_key, n_top_genes, threshold_i, n_perms, n_jobs, **kwargs
        Passed to spatially_variable_genes for each subset.

    Returns
    -------
    dict[str, DataFrame] or None
        Per-library DataFrames (index=genes, columns=spatial_i, spatial_pval, spatially_variable)
        stored in adata.uns['spatial_variable_per_library']. Returns None if library_id_col
        is missing from adata.obs (caller should skip SVG).
    """
    if library_id_col not in adata.obs.columns:
        return None

    libraries = pd.Series(adata.obs[library_id_col]).dropna().unique().tolist()
    if not libraries:
        return None

    out: dict[str, pd.DataFrame] = {}
    for lib in libraries:
        subset = adata[adata.obs[library_id_col] == lib].copy()
        if subset.n_obs < 3 or coord_key not in subset.obsm:
            continue
        try:
            spatially_variable_genes(
                subset,
                mode=mode,
                coord_key=coord_key,
                n_top_genes=n_top_genes,
                threshold_i=threshold_i,
                copy_var=True,
                n_perms=n_perms,
                n_jobs=n_jobs,
                **kwargs,
            )
        except Exception:
            continue
        cols = [c for c in ["spatial_i", "spatial_pval", "spatially_variable"] if c in subset.var.columns]
        if cols:
            out[str(lib)] = subset.var[cols].copy()

    if not out:
        return None
    adata.uns["spatial_variable_per_library"] = out
    return out
