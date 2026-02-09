"""
Signature utility functions.
"""

from typing import List, Optional
import anndata as ad


def get_signature_columns(adata: ad.AnnData, prefix: str = 'sig:', suffix: str = '_z') -> List[str]:
    """
    Get all signature columns from adata.obs.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data object
    prefix : str
        Prefix for signature columns (default: 'sig:')
    suffix : str
        Suffix for signature columns (default: '_z')
    
    Returns
    -------
    list of str
        List of signature column names
    """
    sig_cols = [col for col in adata.obs.columns 
                if col.startswith(prefix) and col.endswith(suffix)]
    return sorted(sig_cols)


def filter_signatures(signature_list: List[str], 
                     include: Optional[List[str]] = None,
                     exclude: Optional[List[str]] = None) -> List[str]:
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
        if normalized.startswith('sig:'):
            normalized = normalized[4:]
        if normalized.endswith('_z'):
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
                if pattern.startswith('sig:'):
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
                if pattern.startswith('sig:'):
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
    cleaned = sig.replace('sig:', '').replace('_z', '')
    # Replace '/' with '-' for readability
    cleaned = cleaned.replace('/', '-')
    # Truncate if too long
    if len(cleaned) > max_length:
        cleaned = cleaned[:max_length-3] + '...'
    return cleaned
