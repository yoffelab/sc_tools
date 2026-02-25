"""
Spatial colocalization analysis utilities.

Provides functions for:
- Pearson correlation
- Truncated similarity (product when both scores > 0)
- Moran's I spatial autocorrelation
- Thresholded neighborhood enrichment
"""

import anndata as ad
import numpy as np
import pandas as pd

from sc_tools.utils.signatures import get_signature_df

try:
    import squidpy as sq

    SQUIDPY_AVAILABLE = True
except ImportError:
    SQUIDPY_AVAILABLE = False
    sq = None

try:
    from tqdm import tqdm
except ImportError:
    # Fallback if tqdm not available
    def tqdm(iterable, desc=None):
        return iterable


def truncated_similarity(score_a, score_b):
    """
    Truncated similarity: score_a * score_b where both > 0, else 0.

    Parameters
    ----------
    score_a : array-like
        First score vector (e.g. proliferation).
    score_b : array-like
        Second score vector (e.g. macrophage).

    Returns
    -------
    np.ndarray
        Same shape as inputs; product where both > 0, else 0.
    """
    a = np.asarray(score_a, dtype=float)
    b = np.asarray(score_b, dtype=float)
    mask = (a > 0) & (b > 0)
    out = np.zeros_like(a)
    out[mask] = a[mask] * b[mask]
    return out


def pearson_correlation(
    adata: ad.AnnData, sig_columns: list[str], min_valid_ratio: float = 0.5
) -> pd.DataFrame:
    """
    Compute Pearson correlation matrix between signatures across spots.

    Parameters
    ----------
    adata : AnnData
        Annotated data object with signature scores
    sig_columns : list of str
        List of signature column names
    min_valid_ratio : float
        Minimum ratio of non-NaN values required (default: 0.5)

    Returns
    -------
    DataFrame
        Correlation matrix
    """
    # Extract signature scores (from obsm or obs)
    sig_df = get_signature_df(adata)
    cols = [c for c in sig_columns if c in sig_df.columns]
    sig_df = sig_df[cols]

    # Remove signatures with too many NaN values
    valid_sigs = sig_df.columns[sig_df.isna().sum() < sig_df.shape[0] * min_valid_ratio].tolist()
    sig_df = sig_df[valid_sigs]

    # Compute correlation matrix
    corr_matrix = sig_df.corr(method="pearson")

    return corr_matrix


def morans_i(
    adata: ad.AnnData,
    sig_column: str,
    coord_key: str = "spatial",
    n_perms: int = 1000,
    n_jobs: int = 1,
) -> dict[str, float]:
    """
    Compute Moran's I spatial autocorrelation for a signature using squidpy.

    Note: squidpy's spatial_autocorr expects genes in var_names, but signature
    scores are typically in obs columns. This function temporarily adds the
    signature as a "pseudo-gene" to compute Moran's I.

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    sig_column : str
        Signature column name in adata.obs
    coord_key : str
        Key in adata.obsm for spatial coordinates (default: 'spatial')
    n_perms : int
        Number of permutations for p-value calculation (default: 1000)
    n_jobs : int
        Number of jobs for parallel processing (default: 1)

    Returns
    -------
    dict
        Dictionary with 'I' (Moran's I) and 'pval' (p-value)
    """
    if not SQUIDPY_AVAILABLE:
        raise ImportError("squidpy is required for Moran's I computation")

    # Ensure spatial neighbors are computed
    if "spatial_neighbors" not in adata.obsp:
        sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)

    # Get signature values
    sig_df = get_signature_df(adata)
    if sig_column not in sig_df.columns and sig_column not in adata.obs.columns:
        raise ValueError(f"Signature column '{sig_column}' not found in obsm or adata.obs")
    sig_values = sig_df[sig_column].values if sig_column in sig_df.columns else adata.obs[sig_column].values

    # Check for NaN values
    if np.isnan(sig_values).all():
        return {"I": np.nan, "pval": np.nan}

    # Create temporary AnnData with signature as a "gene"
    # We need to add it to the expression matrix temporarily
    temp_adata = adata.copy()

    try:
        # Add signature as a temporary "gene" in var
        if sig_column not in temp_adata.var_names:
            # Create a new var entry
            new_var = temp_adata.var.iloc[0:1].copy()
            new_var.index = [sig_column]
            temp_adata.var = pd.concat([temp_adata.var, new_var])

            # Add signature values to X as a new column
            # Convert to dense if sparse
            if hasattr(temp_adata.X, "toarray"):
                X_dense = temp_adata.X.toarray()
            else:
                X_dense = temp_adata.X.copy()

            # Add signature as new column (reshape to column vector)
            sig_col = sig_values.reshape(-1, 1)
            temp_adata.X = np.hstack([X_dense, sig_col])

        # Use squidpy's spatial_autocorr
        sq.gr.spatial_autocorr(
            temp_adata, mode="moran", n_perms=n_perms, n_jobs=n_jobs, genes=[sig_column], copy=False
        )

        # Extract results from uns
        if "moranI" in temp_adata.uns:
            moran_results = temp_adata.uns["moranI"]
            if isinstance(moran_results, pd.DataFrame) and sig_column in moran_results.index:
                I = moran_results.loc[sig_column, "I"]
                pval = moran_results.loc[sig_column, "pval_norm"]
                return {"I": float(I), "pval": float(pval)}

        return {"I": np.nan, "pval": np.nan}

    finally:
        # Restore original state (cleanup)
        # Note: We're working on a copy, so we don't need to restore
        # But it's good practice to be explicit
        pass


def morans_i_batch(
    adata: ad.AnnData,
    sig_columns: list[str],
    n_perms: int = 1000,
    coord_key: str = "spatial",
    n_jobs: int = 1,
) -> pd.DataFrame:
    """
    Compute Moran's I for multiple signatures using squidpy.

    This function efficiently computes Moran's I for multiple signatures by
    temporarily adding them all as "pseudo-genes" and computing in one batch.

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    sig_columns : list of str
        List of signature column names
    n_perms : int
        Number of permutations (default: 1000)
    coord_key : str
        Key in adata.obsm for spatial coordinates (default: 'spatial')
    n_jobs : int
        Number of jobs for parallel processing (default: 1)

    Returns
    -------
    DataFrame
        DataFrame with columns 'Morans_I' and 'p_value'
    """
    if not SQUIDPY_AVAILABLE:
        raise ImportError("squidpy is required for Moran's I computation")

    # Ensure spatial neighbors are computed
    if "spatial_neighbors" not in adata.obsp:
        sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)

    # Filter to valid signatures (non-NaN); get from obsm or obs
    sig_df = get_signature_df(adata)
    valid_sigs = []
    for sig_col in sig_columns:
        if sig_col in sig_df.columns and not sig_df[sig_col].isna().all():
            valid_sigs.append(sig_col)

    if len(valid_sigs) == 0:
        return pd.DataFrame(columns=["Morans_I", "p_value"])

    # Create temporary AnnData with all signatures as "genes"
    temp_adata = adata.copy()

    # Store original state
    original_X = temp_adata.X.copy()
    if hasattr(original_X, "toarray"):
        X_dense = original_X.toarray()
    else:
        X_dense = original_X.copy()

    try:
        # Add all signatures as temporary "genes"
        new_vars = []
        sig_matrix = []

        for sig_col in valid_sigs:
            if sig_col not in temp_adata.var_names:
                # Create var entry
                new_var = temp_adata.var.iloc[0:1].copy()
                new_var.index = [sig_col]
                new_vars.append(new_var)

                # Get signature values (from obsm or obs)
                if sig_col in sig_df.columns:
                    sig_values = sig_df[sig_col].values
                else:
                    sig_values = temp_adata.obs[sig_col].values
                sig_matrix.append(sig_values.reshape(-1, 1))

        if len(new_vars) > 0:
            # Add new var entries
            temp_adata.var = pd.concat([temp_adata.var] + new_vars)

            # Add signature columns to X
            if len(sig_matrix) > 0:
                sig_array = np.hstack(sig_matrix)
                temp_adata.X = np.hstack([X_dense, sig_array])

        # Use squidpy's spatial_autocorr for all signatures at once
        sq.gr.spatial_autocorr(
            temp_adata, mode="moran", n_perms=n_perms, n_jobs=n_jobs, genes=valid_sigs, copy=False
        )

        # Extract results
        results = {}
        if "moranI" in temp_adata.uns:
            moran_results = temp_adata.uns["moranI"]
            if isinstance(moran_results, pd.DataFrame):
                for sig_col in valid_sigs:
                    if sig_col in moran_results.index:
                        I = moran_results.loc[sig_col, "I"]
                        pval = moran_results.loc[sig_col, "pval_norm"]
                        results[sig_col] = {"I": float(I), "pval": float(pval)}
                    else:
                        results[sig_col] = {"I": np.nan, "pval": np.nan}
            else:
                # If results are in different format
                for sig_col in valid_sigs:
                    results[sig_col] = {"I": np.nan, "pval": np.nan}
        else:
            for sig_col in valid_sigs:
                results[sig_col] = {"I": np.nan, "pval": np.nan}

        # Create DataFrame
        moran_df = pd.DataFrame(results).T
        moran_df.columns = ["Morans_I", "p_value"]

        return moran_df

    except Exception as e:
        # Fallback: compute one by one
        print(f"Warning: Batch computation failed, computing individually: {e}")
        results = {}
        for sig_col in tqdm(valid_sigs, desc="Computing Moran's I"):
            try:
                result = morans_i(
                    adata, sig_col, n_perms=n_perms, coord_key=coord_key, n_jobs=n_jobs
                )
                results[sig_col] = result
            except Exception as e2:
                print(f"   Warning: Failed for {sig_col}: {e2}")
                results[sig_col] = {"I": np.nan, "pval": np.nan}

        moran_df = pd.DataFrame(results).T
        moran_df.columns = ["Morans_I", "p_value"]
        return moran_df


def neighborhood_enrichment(
    adata: ad.AnnData,
    sig_column: str,
    threshold_low: float = -1.0,
    threshold_high: float = 1.0,
    n_perms: int = 1000,
) -> float:
    """
    Compute thresholded neighborhood enrichment for a signature.

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    sig_column : str
        Signature column name
    threshold_low : float
        Threshold for low category (default: -1.0)
    threshold_high : float
        Threshold for high category (default: 1.0)
    n_perms : int
        Number of permutations (default: 1000)

    Returns
    -------
    float
        Neighborhood enrichment score
    """
    # Threshold signature (from obsm or obs)
    sig_df = get_signature_df(adata)
    if sig_column in sig_df.columns:
        scores = sig_df[sig_column].values
    elif sig_column in adata.obs.columns:
        scores = adata.obs[sig_column].values
    else:
        raise ValueError(f"Signature column '{sig_column}' not found in obsm or adata.obs")
    cat_binary = (scores >= threshold_high).astype(str)

    # Store as categorical
    temp_key = f"_temp_{sig_column}"
    adata.obs[temp_key] = pd.Categorical(cat_binary)

    # Ensure spatial neighbors are computed
    if "spatial_neighbors" not in adata.obsp:
        if not SQUIDPY_AVAILABLE:
            raise ImportError("squidpy is required for neighborhood enrichment")
        sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)

    if not SQUIDPY_AVAILABLE:
        raise ImportError("squidpy is required for neighborhood enrichment")

    try:
        # Compute neighborhood enrichment
        sq.gr.nhood_enrichment(
            adata,
            cluster_key=temp_key,
            n_perms=n_perms,
        )

        if "nhood_enrichment" in adata.uns:
            enrichment = adata.uns["nhood_enrichment"]
            if (
                isinstance(enrichment, pd.DataFrame)
                and "True" in enrichment.index
                and "True" in enrichment.columns
            ):
                result = enrichment.loc["True", "True"]
            else:
                result = np.nan
        else:
            result = np.nan
    except Exception:
        result = np.nan
    finally:
        # Clean up
        if temp_key in adata.obs.columns:
            del adata.obs[temp_key]

    return result


def neighborhood_enrichment_batch(
    adata: ad.AnnData,
    sig_columns: list[str],
    threshold_low: float = -1.0,
    threshold_high: float = 1.0,
    n_perms: int = 1000,
) -> pd.DataFrame:
    """
    Compute neighborhood enrichment for multiple signatures.

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    sig_columns : list of str
        List of signature column names
    threshold_low : float
        Threshold for low category (default: -1.0)
    threshold_high : float
        Threshold for high category (default: 1.0)
    n_perms : int
        Number of permutations (default: 1000)

    Returns
    -------
    DataFrame
        DataFrame with 'nhood_enrichment' column
    """
    results = {}

    for sig in tqdm(sig_columns, desc="Computing neighborhood enrichment"):
        try:
            enrichment = neighborhood_enrichment(
                adata,
                sig,
                threshold_low=threshold_low,
                threshold_high=threshold_high,
                n_perms=n_perms,
            )
            results[sig] = enrichment
        except Exception as e:
            print(f"   Warning: Neighborhood enrichment failed for {sig}: {e}")
            results[sig] = np.nan

    # Create DataFrame
    nhood_df = pd.DataFrame.from_dict(results, orient="index", columns=["nhood_enrichment"])

    return nhood_df


__all__ = [
    "truncated_similarity",
    "pearson_correlation",
    "morans_i",
    "morans_i_batch",
    "neighborhood_enrichment",
    "neighborhood_enrichment_batch",
]
