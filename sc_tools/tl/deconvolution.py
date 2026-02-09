"""
Deconvolution utilities: signature gene selection and caching.
"""

import logging
from pathlib import Path
from typing import Optional, List
import scanpy as sc
import numpy as np
import anndata as ad

from ..data.io import get_cache_key, load_cached_signatures, save_cached_signatures

logger = logging.getLogger(__name__)


def select_signature_genes(sc_adata: ad.AnnData, celltype_key: str, 
                           sc_batch_key: str, n_genes_max: int, 
                           skip_hvg: bool = True, 
                           cache_dir: Optional[str] = None,
                           force_recompute: bool = False,
                           sc_data_file: Optional[str] = None,
                           logger_instance: Optional[logging.Logger] = None) -> List[str]:
    """
    Select signature genes for deconvolution.
    
    This function selects genes for deconvolution by combining:
    1. Highly variable genes (HVGs) - optional
    2. Top marker genes per cell type (always computed)
    
    Results can be cached to avoid recomputation.
    
    Parameters
    ----------
    sc_adata : AnnData
        Single-cell reference data
    celltype_key : str
        Column name in sc_adata.obs containing cell type annotations
    sc_batch_key : str
        Column name in sc_adata.obs containing batch information
    n_genes_max : int
        Maximum number of genes to return
    skip_hvg : bool
        If True, skip HVG computation and use only marker genes (faster, less memory)
    cache_dir : str, optional
        Directory to cache signature genes. If None, caching is disabled.
    force_recompute : bool
        If True, force recomputation even if cache exists.
    sc_data_file : str, optional
        Path to single-cell data file (for cache key generation).
    logger_instance : Logger, optional
        Custom logger instance. If None, uses module logger.
    
    Returns
    -------
    list of str
        List of signature gene names
    """
    log = logger_instance if logger_instance is not None else logger
    log.info("Selecting signature genes...")
    
    # Check cache if enabled
    if cache_dir and not force_recompute and sc_data_file:
        cache_key = get_cache_key(sc_data_file, celltype_key, sc_batch_key, n_genes_max, skip_hvg)
        cache_path = Path(cache_dir) / f"signature_genes_{cache_key}.pkl"
        
        cached_genes = load_cached_signatures(cache_path, logger_instance=log)
        if cached_genes is not None:
            return cached_genes
        log.info("   Cache miss or invalid, computing signature genes...")
    elif force_recompute:
        log.info("   Force recompute enabled, computing signature genes...")
    
    candidate_genes = set()
    
    # Step 1: Check if HVGs are already computed
    if not skip_hvg and 'highly_variable' in sc_adata.var.columns:
        log.info("   Using pre-computed highly variable genes...")
        hvg_genes = set(sc_adata.var[sc_adata.var.highly_variable].index)
        candidate_genes.update(hvg_genes)
        log.info(f"   Found {len(hvg_genes)} pre-computed HVGs")
    elif not skip_hvg:
        # Try to compute HVGs with error handling
        try:
            log.info("   Computing highly variable genes...")
            # Use a simpler method that's less memory intensive
            sc.pp.highly_variable_genes(sc_adata, flavor="seurat", n_top_genes=2000, 
                                       subset=False, batch_key=sc_batch_key)
            hvg_genes = set(sc_adata.var[sc_adata.var.highly_variable].index)
            candidate_genes.update(hvg_genes)
            log.info(f"   Found {len(hvg_genes)} HVGs")
        except Exception as e:
            log.warning(f"   Failed to compute HVGs: {e}")
            log.warning("   Continuing with marker genes only...")
    
    # Step 2: Top marker genes per cell type (always compute these)
    log.info("   Computing marker genes per cell type...")
    try:
        sc.tl.rank_genes_groups(sc_adata, groupby=celltype_key, method='wilcoxon', use_raw=False)
        marker_genes = []
        for group in sc_adata.obs[celltype_key].unique():
            df = sc.get.rank_genes_groups_df(sc_adata, group=group)
            top_genes = df.head(100).names.tolist()
            marker_genes.extend(top_genes)
        marker_genes = set(marker_genes)
        candidate_genes.update(marker_genes)
        log.info(f"   Found {len(marker_genes)} marker genes")
    except Exception as e:
        log.error(f"   Failed to compute marker genes: {e}")
        raise
    
    # Step 3: If no genes selected, use top expressed genes as fallback
    if len(candidate_genes) == 0:
        log.warning("   No genes selected, using top expressed genes as fallback...")
        # Calculate mean expression per gene
        if hasattr(sc_adata.X, 'toarray'):
            mean_expr = np.array(sc_adata.X.mean(axis=0)).flatten()
        else:
            mean_expr = np.array(sc_adata.X.mean(axis=0)).flatten()
        top_indices = np.argsort(mean_expr)[-n_genes_max:][::-1]
        candidate_genes = set(sc_adata.var_names[top_indices])
        log.info(f"   Selected {len(candidate_genes)} top expressed genes")
    
    # Step 4: Limit to max genes
    candidate_genes = list(candidate_genes)
    if len(candidate_genes) > n_genes_max:
        # Prioritize marker genes if we have both
        if 'marker_genes' in locals() and len(marker_genes) > 0:
            # Keep all marker genes, then fill with HVGs or top expressed
            marker_list = list(marker_genes)
            other_genes = [g for g in candidate_genes if g not in marker_genes]
            candidate_genes = marker_list + other_genes[:n_genes_max - len(marker_list)]
        else:
            candidate_genes = candidate_genes[:n_genes_max]
    
    log.info(f"   Using {len(candidate_genes)} signature genes")
    
    # Save to cache if enabled
    if cache_dir and sc_data_file:
        cache_key = get_cache_key(sc_data_file, celltype_key, sc_batch_key, n_genes_max, skip_hvg)
        cache_path = Path(cache_dir) / f"signature_genes_{cache_key}.pkl"
        metadata = {
            'sc_data_file': sc_data_file,
            'celltype_key': celltype_key,
            'sc_batch_key': sc_batch_key,
            'n_genes_max': n_genes_max,
            'skip_hvg': skip_hvg,
            'n_genes': len(candidate_genes)
        }
        save_cached_signatures(candidate_genes, cache_path, metadata, logger_instance=log)
    
    return candidate_genes
