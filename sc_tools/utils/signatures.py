"""
Signature utility functions.
"""

from __future__ import annotations

import anndata as ad
import pandas as pd


def get_signature_columns_from_obsm(
    adata: ad.AnnData,
    obsm_key: str = "signature_score_z",
) -> list[str]:
    """
    Get signature column names from adata.obsm (full-path style, e.g. Myeloid/Macrophage_Core).

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    obsm_key : str
        Key in adata.obsm (default: 'signature_score_z')

    Returns
    -------
    list of str
        Column names, or empty list if key is missing
    """
    if obsm_key not in adata.obsm:
        return []
    arr = adata.obsm[obsm_key]
    if hasattr(arr, "columns"):
        return list(arr.columns)
    return []


def get_signature_columns(adata: ad.AnnData, prefix: str = "sig:", suffix: str = "_z") -> list[str]:
    """
    Get all signature column names. Prefer obsm['signature_score_z'] if present; else obs columns.

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    prefix : str
        Prefix for obs-based signature columns (default: 'sig:')
    suffix : str
        Suffix for obs-based signature columns (default: '_z')

    Returns
    -------
    list of str
        List of signature column names (obsm full-path or obs sig:..._z)
    """
    obsm_cols = get_signature_columns_from_obsm(adata, obsm_key="signature_score_z")
    if obsm_cols:
        return sorted(obsm_cols)
    sig_cols = [col for col in adata.obs.columns if col.startswith(prefix) and col.endswith(suffix)]
    return sorted(sig_cols)


def get_signature_df(
    adata: ad.AnnData,
    use_z: bool = True,
    obsm_key: str | None = None,
) -> pd.DataFrame:
    """
    Return a DataFrame of signature scores. Prefer obsm; fall back to obs for backward compatibility.

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    use_z : bool
        If True (default), use z-scored scores (signature_score_z); else signature_score
    obsm_key : str or None
        If set, use this obsm key; else use signature_score_z or signature_score

    Returns
    -------
    DataFrame
        Index = adata.obs_names; columns = signature names (full-path or obs column names)
    """
    key = obsm_key
    if key is None:
        key = "signature_score_z" if use_z else "signature_score"
    if key in adata.obsm:
        df = adata.obsm[key]
        if hasattr(df, "copy"):
            return df.copy()
        return pd.DataFrame(df, index=adata.obs_names)
    # Fallback: build from obs (prefix sig:, suffix _z)
    prefix, suffix = "sig:", "_z"
    cols = [c for c in adata.obs.columns if c.startswith(prefix) and c.endswith(suffix)]
    if not cols:
        return pd.DataFrame(index=adata.obs_names)
    return adata.obs[cols].copy()


def filter_signatures(
    signature_list: list[str],
    include: list[str] | None = None,
    exclude: list[str] | None = None,
) -> list[str]:
    """
    Filter signatures based on include/exclude criteria.

    Parameters
    ----------
    signature_list : list
        List of signature column names (e.g., 'sig:Tumor_Cells-EMT_Tumor_z')
    include : list or None
        List of signatures to include. Can be:
        - Exact matches (with or without 'sig:' prefix and '_z' suffix)
        - Patterns (substring search, case-insensitive)
        - None to include all
    exclude : list or None
        List of signatures to exclude. Same format as include.
        - None to exclude none

    Returns
    -------
    list
        Filtered list of signatures
    """
    if include is None and exclude is None:
        return signature_list

    filtered = []

    # Normalize signature names for matching (remove prefix/suffix)
    def normalize_sig(sig):
        """Remove prefix and suffix for comparison."""
        normalized = sig
        if normalized.startswith("sig:"):
            normalized = normalized[4:]
        if normalized.endswith("_z"):
            normalized = normalized[:-2]
        return normalized.lower()

    # Normalize all signatures
    sig_normalized = {sig: normalize_sig(sig) for sig in signature_list}

    for sig in signature_list:
        sig_norm = sig_normalized[sig]
        include_match = False
        exclude_match = False

        # Check include list
        if include is not None:
            for pattern in include:
                # Normalize pattern
                if pattern.startswith("sig:"):
                    pattern_norm = normalize_sig(pattern)
                else:
                    pattern_norm = pattern.lower()

                # Check exact match (normalized) or substring match
                if pattern_norm == sig_norm or pattern_norm in sig_norm:
                    include_match = True
                    break
        else:
            # No include list means include all
            include_match = True

        # Check exclude list
        if exclude is not None:
            for pattern in exclude:
                # Normalize pattern
                if pattern.startswith("sig:"):
                    pattern_norm = normalize_sig(pattern)
                else:
                    pattern_norm = pattern.lower()

                # Check exact match (normalized) or substring match
                if pattern_norm == sig_norm or pattern_norm in sig_norm:
                    exclude_match = True
                    break

        # Include if matches include criteria and doesn't match exclude
        if include_match and not exclude_match:
            filtered.append(sig)

    return filtered


def clean_sig_name(sig: str, max_length: int = 40) -> str:
    """
    Clean signature name for display.

    Parameters
    ----------
    sig : str
        Signature name (e.g., 'sig:Tumor_Cells/EMT_Tumor_z')
    max_length : int
        Maximum length for truncated names (default: 40)

    Returns
    -------
    str
        Cleaned signature name
    """
    # Remove 'sig:' prefix and '_z' suffix
    cleaned = sig.replace("sig:", "").replace("_z", "")
    # Replace '/' with '-' for readability
    cleaned = cleaned.replace("/", "-")
    # Truncate if too long
    if len(cleaned) > max_length:
        cleaned = cleaned[: max_length - 3] + "..."
    return cleaned
