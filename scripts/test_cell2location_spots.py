"""
Test Cell2location with progressively larger spot counts.

This script tests Cell2location deconvolution with different numbers of spots
(1, 10, 100, 1000) to identify when the program fails or encounters issues.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import os
import gc
import pickle
import hashlib
from pathlib import Path
from tqdm import tqdm
from typing import Optional, List, Dict
import warnings
import sys
import logging
from datetime import datetime
import traceback
import signal
import time

# Allow mutex warnings to be visible for debugging
# os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'  # Commented out to see all warnings
os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'max_split_size_mb:128'  # Reduce CUDA memory fragmentation

# Note: Mutex warnings are now visible. They come from C++ code in TensorFlow/PyTorch
# and indicate thread synchronization, which is usually harmless but can indicate
# performance issues or deadlocks if excessive.

# Memory profiling
try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False

# Keep warnings visible for debugging (can re-enable filtering if needed)
# warnings.filterwarnings('ignore', category=UserWarning)
# warnings.filterwarnings('ignore', category=FutureWarning)
# warnings.filterwarnings('ignore', category=DeprecationWarning)

# Allow TensorFlow/PyTorch verbose output to see mutex warnings
# Commented out to see all warnings including mutex warnings
# try:
#     import tensorflow as tf
#     tf.get_logger().setLevel('ERROR')
# except ImportError:
#     pass

# try:
#     import torch
#     import logging as torch_logging
#     torch_logging.getLogger('torch').setLevel('ERROR')
# except ImportError:
#     pass

# Setup logging
log_dir = 'output/deconvolution/test_logs'
os.makedirs(log_dir, exist_ok=True)
log_file = os.path.join(log_dir, f'cell2location_test_{datetime.now().strftime("%Y%m%d_%H%M%S")}.txt')
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

sc.settings.set_figure_params(dpi=300, dpi_save=400)

# ============================================================================
# Configuration
# ============================================================================
CONFIG = {
    'sc_data_file': "results/seurat_object.h5ad",
    'spatial_data_file': "results/adata.img.genescores.h5ad",
    'celltype_key': 'cell.type',
    'spatial_batch_key': 'library_id',
    'sc_batch_key': 'Batch',
    'n_genes_max': 2000,
    'cell2location_max_epochs': 100,  # Very reduced for testing (was 25000)
    'test_spot_counts': [1, 10, 100, 1000],  # Number of spots to test
    'skip_hvg_computation': True,  # Skip HVG computation to avoid memory issues
    'max_time_per_test_minutes': 30,  # Maximum time per test in minutes (safety timeout)
    'use_gpu': None,  # None = auto-detect GPU, True = force GPU, False = force CPU
    'cache_signature_genes': True,  # Cache signature genes to avoid recomputation
    'signature_cache_dir': 'output/deconvolution/cache',  # Directory for cache files
    'force_recompute_signatures': False,  # Force recomputation even if cache exists
}

# ============================================================================
# Helper Functions
# ============================================================================

def get_memory_usage():
    """Get current memory usage in MB."""
    if PSUTIL_AVAILABLE:
        process = psutil.Process(os.getpid())
        mem_info = process.memory_info()
        return mem_info.rss / 1024 / 1024  # RSS in MB
    return 0.0


def log_memory(step_name: str):
    """Log memory usage at a specific step."""
    mem_mb = get_memory_usage()
    logger.info(f"[MEMORY] {step_name}: RSS={mem_mb:.1f}MB")


def aggressive_cleanup():
    """Aggressively clean up memory."""
    gc.collect()
    gc.collect()  # Call twice
    if hasattr(gc, 'collect'):
        try:
            import ctypes
            libc = ctypes.CDLL("libc.dylib" if sys.platform == "darwin" else "libc.so.6")
            libc.malloc_trim(0)
        except:
            pass


def get_cache_key(sc_data_file: str, celltype_key: str, sc_batch_key: str, 
                 n_genes_max: int, skip_hvg: bool) -> str:
    """Generate a cache key based on input parameters and file modification time."""
    # Include file modification time to detect changes
    try:
        mtime = os.path.getmtime(sc_data_file)
    except OSError:
        mtime = 0
    
    # Create a hash from parameters
    key_string = f"{sc_data_file}:{mtime}:{celltype_key}:{sc_batch_key}:{n_genes_max}:{skip_hvg}"
    key_hash = hashlib.md5(key_string.encode()).hexdigest()
    return key_hash


def load_cached_signatures(cache_path: Path) -> Optional[List[str]]:
    """Load cached signature genes from file."""
    try:
        if cache_path.exists():
            with open(cache_path, 'rb') as f:
                cache_data = pickle.load(f)
                if isinstance(cache_data, dict) and 'signature_genes' in cache_data:
                    logger.info(f"   ✅ Loaded {len(cache_data['signature_genes'])} signature genes from cache")
                    return cache_data['signature_genes']
                elif isinstance(cache_data, list):
                    # Backward compatibility: old format was just a list
                    logger.info(f"   ✅ Loaded {len(cache_data)} signature genes from cache (old format)")
                    return cache_data
    except Exception as e:
        logger.warning(f"   ⚠️  Failed to load cache: {e}")
    return None


def save_cached_signatures(signature_genes: List[str], cache_path: Path, 
                           metadata: Dict) -> None:
    """Save signature genes to cache file with metadata."""
    try:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_data = {
            'signature_genes': signature_genes,
            'metadata': metadata,
            'timestamp': datetime.now().isoformat()
        }
        with open(cache_path, 'wb') as f:
            pickle.dump(cache_data, f)
        logger.info(f"   💾 Saved {len(signature_genes)} signature genes to cache: {cache_path}")
    except Exception as e:
        logger.warning(f"   ⚠️  Failed to save cache: {e}")


def load_single_cell_data(sc_data_file: str, celltype_key: str, 
                          sc_batch_key: str, qc_labels: list) -> tuple:
    """Load and preprocess single-cell reference data."""
    logger.info("="*70)
    logger.info("Loading single-cell reference data")
    logger.info("="*70)
    
    logger.info(f"Loading from: {sc_data_file}")
    sc_adata = sc.read_h5ad(sc_data_file)
    sc_adata.var_names_make_unique()
    log_memory("After loading single-cell data")
    logger.info(f"Single-cell data: {sc_adata.shape} (cells x genes)")
    
    # Preprocess
    if np.max(sc_adata.X) > 100:  # Likely raw counts
        logger.info("Normalizing and log-transforming scRNA-seq data...")
        sc.pp.normalize_total(sc_adata, target_sum=1e4)
        sc.pp.log1p(sc_adata)
    else:
        logger.info("Single-cell data already appears log-normalized")
    
    # Get cell type names (exclude QC labels)
    if celltype_key in sc_adata.obs.columns:
        unique_cts = sorted(sc_adata.obs[celltype_key].unique())
        unique_cts = [ct for ct in unique_cts if ct not in qc_labels]
        logger.info(f"Cell types (excluding QC labels): {len(unique_cts)}")
    else:
        raise ValueError(f"celltype_key '{celltype_key}' not found in single-cell data")
    
    return sc_adata, unique_cts


def select_signature_genes(sc_adata: ad.AnnData, celltype_key: str, 
                           sc_batch_key: str, n_genes_max: int, 
                           skip_hvg: bool = True, 
                           cache_dir: Optional[str] = None,
                           force_recompute: bool = False,
                           sc_data_file: Optional[str] = None) -> list:
    """
    Select signature genes for deconvolution.
    
    Parameters
    ----------
    skip_hvg : bool
        If True, skip HVG computation and use only marker genes (faster, less memory)
    cache_dir : str, optional
        Directory to cache signature genes. If None, caching is disabled.
    force_recompute : bool
        If True, force recomputation even if cache exists.
    sc_data_file : str, optional
        Path to single-cell data file (for cache key generation).
    """
    logger.info("Selecting signature genes...")
    
    # Check cache if enabled
    if cache_dir and not force_recompute and sc_data_file:
        cache_key = get_cache_key(sc_data_file, celltype_key, sc_batch_key, n_genes_max, skip_hvg)
        cache_path = Path(cache_dir) / f"signature_genes_{cache_key}.pkl"
        
        cached_genes = load_cached_signatures(cache_path)
        if cached_genes is not None:
            return cached_genes
        logger.info("   Cache miss or invalid, computing signature genes...")
    elif force_recompute:
        logger.info("   Force recompute enabled, computing signature genes...")
    
    candidate_genes = set()
    
    # Step 1: Check if HVGs are already computed
    if not skip_hvg and 'highly_variable' in sc_adata.var.columns:
        logger.info("   Using pre-computed highly variable genes...")
        hvg_genes = set(sc_adata.var[sc_adata.var.highly_variable].index)
        candidate_genes.update(hvg_genes)
        logger.info(f"   Found {len(hvg_genes)} pre-computed HVGs")
    elif not skip_hvg:
        # Try to compute HVGs with error handling
        try:
            logger.info("   Computing highly variable genes...")
            # Use a simpler method that's less memory intensive
            sc.pp.highly_variable_genes(sc_adata, flavor="seurat", n_top_genes=2000, 
                                       subset=False, batch_key=sc_batch_key)
            hvg_genes = set(sc_adata.var[sc_adata.var.highly_variable].index)
            candidate_genes.update(hvg_genes)
            logger.info(f"   Found {len(hvg_genes)} HVGs")
        except Exception as e:
            logger.warning(f"   Failed to compute HVGs: {e}")
            logger.warning("   Continuing with marker genes only...")
    
    # Step 2: Top marker genes per cell type (always compute these)
    logger.info("   Computing marker genes per cell type...")
    try:
        sc.tl.rank_genes_groups(sc_adata, groupby=celltype_key, method='wilcoxon', use_raw=False)
        marker_genes = []
        for group in sc_adata.obs[celltype_key].unique():
            df = sc.get.rank_genes_groups_df(sc_adata, group=group)
            top_genes = df.head(100).names.tolist()
            marker_genes.extend(top_genes)
        marker_genes = set(marker_genes)
        candidate_genes.update(marker_genes)
        logger.info(f"   Found {len(marker_genes)} marker genes")
    except Exception as e:
        logger.error(f"   Failed to compute marker genes: {e}")
        raise
    
    # Step 3: If no genes selected, use top expressed genes as fallback
    if len(candidate_genes) == 0:
        logger.warning("   No genes selected, using top expressed genes as fallback...")
        # Calculate mean expression per gene
        if hasattr(sc_adata.X, 'toarray'):
            mean_expr = np.array(sc_adata.X.mean(axis=0)).flatten()
        else:
            mean_expr = np.array(sc_adata.X.mean(axis=0)).flatten()
        top_indices = np.argsort(mean_expr)[-n_genes_max:][::-1]
        candidate_genes = set(sc_adata.var_names[top_indices])
        logger.info(f"   Selected {len(candidate_genes)} top expressed genes")
    
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
    
    logger.info(f"   Using {len(candidate_genes)} signature genes")
    
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
        save_cached_signatures(candidate_genes, cache_path, metadata)
    
    return candidate_genes


class TimeoutError(Exception):
    """Custom timeout exception."""
    pass

def timeout_handler(signum, frame):
    """Handler for timeout signal."""
    raise TimeoutError("Operation timed out")

def check_gpu_available() -> bool:
    """Check if GPU is available for PyTorch."""
    try:
        import torch
        if torch.cuda.is_available():
            logger.info(f"   ✅ GPU available: {torch.cuda.get_device_name(0)}")
            logger.info(f"   GPU memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.2f} GB")
            return True
        else:
            logger.info("   ⚠️  GPU not available, using CPU")
            return False
    except ImportError:
        logger.info("   ⚠️  PyTorch not available, using CPU")
        return False
    except Exception as e:
        logger.warning(f"   ⚠️  Error checking GPU: {e}, using CPU")
        return False

def test_cell2location(sc_adata: ad.AnnData, spatial_subset: ad.AnnData,
                       shared_genes: list, celltype_key: str,
                       max_epochs: int, n_spots: int, 
                       max_time_minutes: int = 30, use_gpu: bool = None) -> dict:
    """
    Test Cell2location deconvolution on a spatial subset.
    
    Returns
    -------
    dict
        Results dictionary with 'success', 'error', 'proportions', 'memory_peak'
    """
    result = {
        'success': False,
        'error': None,
        'error_type': None,
        'proportions_shape': None,
        'memory_peak_mb': 0.0,
        'memory_final_mb': 0.0,
    }
    
    try:
        import cell2location as c2l
        from cell2location.models import RegressionModel, Cell2location
        
        logger.info(f"Testing Cell2location with {n_spots} spots...")
        logger.info(f"Maximum time allowed: {max_time_minutes} minutes")
        
        # Check GPU availability
        if use_gpu is None:
            use_gpu = check_gpu_available()
        elif use_gpu:
            # User explicitly requested GPU, verify it's available
            if not check_gpu_available():
                logger.warning("   GPU requested but not available, falling back to CPU")
                use_gpu = False
        else:
            logger.info("   Using CPU (use_gpu=False)")
        
        start_time = time.time()
        log_memory(f"Start ({n_spots} spots)")
        
        # Prepare data - Cell2location needs raw counts
        logger.info("   Creating gene subsets...")
        sc_copy = sc_adata[:, shared_genes].copy()
        spatial_subset_sig = spatial_subset[:, shared_genes].copy()
        log_memory(f"After creating subsets ({n_spots} spots)")
        
        # Ensure raw counts
        if sc_copy.raw is None:
            sc_copy.raw = sc_copy.copy()
        if spatial_subset_sig.raw is None:
            spatial_subset_sig.raw = spatial_subset_sig.copy()
        
        # Train regression model on single-cell data
        logger.info("   Training regression model on single-cell data...")
        c2l.models.RegressionModel.setup_anndata(
            sc_copy, labels_key=celltype_key, batch_key=None
        )
        mod = RegressionModel(sc_copy)
        log_memory(f"After creating regression model ({n_spots} spots)")
        
        # Use very reduced parameters for faster testing
        actual_epochs = min(max_epochs, 100)  # Very reduced for testing (was 1000)
        num_samples = 50  # Reduced from 100
        batch_size = 2500
        
        logger.info(f"   Training regression model (epochs={actual_epochs})...")
        logger.info(f"   This may take a few minutes...")
        
        # Add progress callback if available
        try:
            mod.train(max_epochs=actual_epochs, use_gpu=use_gpu, 
                     progress_bar_refresh_rate=10)  # Show progress every 10 epochs
        except TypeError:
            # Fallback if progress_bar_refresh_rate not supported
            mod.train(max_epochs=actual_epochs, use_gpu=use_gpu)
        
        log_memory(f"After training regression model ({n_spots} spots)")
        
        logger.info("   Exporting posterior from single-cell model...")
        sc_copy = mod.export_posterior(sc_copy, sample_kwargs={'num_samples': num_samples, 'batch_size': batch_size})
        log_memory(f"After exporting posterior ({n_spots} spots)")
        
        # Train Cell2location on spatial data
        logger.info("   Training on spatial data...")
        c2l.models.Cell2location.setup_anndata(
            spatial_subset_sig, batch_key=None
        )
        mod_spatial = Cell2location(spatial_subset_sig, sc_copy)
        log_memory(f"After creating spatial model ({n_spots} spots)")
        
        logger.info(f"   Training spatial model (epochs={actual_epochs})...")
        logger.info(f"   This may take a few minutes...")
        
        # Add progress callback if available
        try:
            mod_spatial.train(max_epochs=actual_epochs, use_gpu=use_gpu,
                             progress_bar_refresh_rate=10)  # Show progress every 10 epochs
        except TypeError:
            # Fallback if progress_bar_refresh_rate not supported
            mod_spatial.train(max_epochs=actual_epochs, use_gpu=use_gpu)
        
        log_memory(f"After training spatial model ({n_spots} spots)")
        
        # Get proportions
        logger.info("   Extracting proportions...")
        spatial_subset_sig = mod_spatial.export_posterior(
            spatial_subset_sig, sample_kwargs={'num_samples': num_samples, 'batch_size': batch_size}
        )
        log_memory(f"After exporting spatial posterior ({n_spots} spots)")
        
        # Extract proportions from obsm
        if 'q05_cell_abundance_w_sf' in spatial_subset_sig.obsm:
            proportions = spatial_subset_sig.obsm['q05_cell_abundance_w_sf'].values.copy()
        elif 'means_cell_abundance_w_sf' in spatial_subset_sig.obsm:
            proportions = spatial_subset_sig.obsm['means_cell_abundance_w_sf'].values.copy()
        else:
            raise ValueError(f"Could not find Cell2location proportions in obsm. Available keys: {list(spatial_subset_sig.obsm.keys())}")
        
        # Normalize to proportions
        row_sums = proportions.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1
        proportions = proportions / row_sums
        
        result['success'] = True
        result['proportions_shape'] = proportions.shape
        result['memory_final_mb'] = get_memory_usage()
        
        # Cleanup
        mod = None
        mod_spatial = None
        sc_copy = None
        spatial_subset_sig = None
        aggressive_cleanup()
        
        logger.info(f"✅ SUCCESS with {n_spots} spots! Proportions shape: {proportions.shape}")
        return result
            
    except ImportError:
        result['error'] = "Cell2location not installed"
        result['error_type'] = "ImportError"
        logger.error("❌ Cell2location not installed")
        return result
    except Exception as e:
        result['error'] = str(e)
        result['error_type'] = type(e).__name__
        logger.error(f"❌ ERROR with {n_spots} spots: {e}")
        logger.error(traceback.format_exc())
        aggressive_cleanup()
        return result


def main():
    """Main function to test Cell2location with different spot counts."""
    logger.info("="*70)
    logger.info("Cell2location Spot Count Testing")
    logger.info("="*70)
    
    # Load single-cell data
    logger.info("\n1. Loading single-cell reference data...")
    qc_labels = ['QC', 'Doublet', 'Ambiguous']
    sc_adata, unique_cts = load_single_cell_data(
        CONFIG['sc_data_file'],
        CONFIG['celltype_key'],
        CONFIG['sc_batch_key'],
        qc_labels
    )
    
    # Select signature genes (skip HVG computation to avoid memory issues)
    logger.info("\n2. Selecting signature genes...")
    logger.info("   (Skipping HVG computation to avoid memory issues)")
    shared_genes = select_signature_genes(
        sc_adata,
        CONFIG['celltype_key'],
        CONFIG['sc_batch_key'],
        CONFIG['n_genes_max'],
        skip_hvg=CONFIG.get('skip_hvg_computation', True),  # Skip HVG computation
        cache_dir=CONFIG.get('signature_cache_dir') if CONFIG.get('cache_signature_genes', False) else None,
        force_recompute=CONFIG.get('force_recompute_signatures', False),
        sc_data_file=CONFIG['sc_data_file']
    )
    
    # Ensure shared genes are in both datasets
    sc_adata = sc_adata[:, shared_genes].copy()
    
    # Load spatial data
    logger.info("\n3. Loading spatial data...")
    spatial_adata = sc.read_h5ad(CONFIG['spatial_data_file'])
    spatial_adata.var_names_make_unique()
    
    # Get first library for testing
    library_ids = sorted(spatial_adata.obs[CONFIG['spatial_batch_key']].unique())
    if len(library_ids) == 0:
        raise ValueError("No libraries found in spatial data")
    
    test_library = library_ids[0]
    logger.info(f"Using library: {test_library}")
    
    spatial_lib = spatial_adata[spatial_adata.obs[CONFIG['spatial_batch_key']] == test_library].copy()
    logger.info(f"Spatial data shape: {spatial_lib.shape}")
    
    # Find shared genes
    shared_genes_spatial = [g for g in shared_genes if g in spatial_lib.var_names]
    logger.info(f"Shared genes: {len(shared_genes_spatial)}")
    
    if len(shared_genes_spatial) < 100:
        raise ValueError(f"Too few shared genes: {len(shared_genes_spatial)}")
    
    # Test with different spot counts
    logger.info("\n" + "="*70)
    logger.info("Testing Cell2location with different spot counts")
    logger.info("="*70)
    
    results = []
    
    for n_spots in CONFIG['test_spot_counts']:
        logger.info("\n" + "-"*70)
        logger.info(f"TESTING WITH {n_spots} SPOTS")
        logger.info("-"*70)
        
        # Subset spatial data
        if n_spots > len(spatial_lib):
            logger.warning(f"Requested {n_spots} spots but only {len(spatial_lib)} available. Using all spots.")
            n_spots = len(spatial_lib)
        
        # Randomly sample spots
        np.random.seed(42)  # For reproducibility
        spot_indices = np.random.choice(len(spatial_lib), size=n_spots, replace=False)
        spatial_subset = spatial_lib[spot_indices].copy()
        
        logger.info(f"Subset shape: {spatial_subset.shape} (spots x genes)")
        log_memory(f"After subsetting to {n_spots} spots")
        
        # Test Cell2location
        result = test_cell2location(
            sc_adata,
            spatial_subset,
            shared_genes_spatial,
            CONFIG['celltype_key'],
            CONFIG['cell2location_max_epochs'],
            n_spots,
            max_time_minutes=CONFIG.get('max_time_per_test_minutes', 30),
            use_gpu=CONFIG.get('use_gpu', None)  # None = auto-detect
        )
        
        result['n_spots'] = n_spots
        result['n_genes'] = len(shared_genes_spatial)
        results.append(result)
        
        # Cleanup
        spatial_subset = None
        aggressive_cleanup()
        log_memory(f"After cleanup ({n_spots} spots)")
    
    # Summary
    logger.info("\n" + "="*70)
    logger.info("TEST SUMMARY")
    logger.info("="*70)
    
    summary_df = pd.DataFrame(results)
    print("\nResults:")
    print(summary_df[['n_spots', 'success', 'error_type', 'proportions_shape', 'memory_final_mb']].to_string())
    
    # Save results
    output_csv = os.path.join(log_dir, 'cell2location_test_results.csv')
    summary_df.to_csv(output_csv, index=False)
    logger.info(f"\nResults saved to: {output_csv}")
    
    # Check which tests passed
    successful = summary_df[summary_df['success'] == True]
    failed = summary_df[summary_df['success'] == False]
    
    logger.info(f"\n✅ Successful tests: {len(successful)}/{len(results)}")
    if len(successful) > 0:
        logger.info(f"   Largest successful: {successful['n_spots'].max()} spots")
    
    if len(failed) > 0:
        logger.info(f"\n❌ Failed tests: {len(failed)}/{len(results)}")
        for _, row in failed.iterrows():
            logger.info(f"   {row['n_spots']} spots: {row['error_type']} - {row['error']}")
    
    logger.info(f"\nFull log saved to: {log_file}")


if __name__ == "__main__":
    main()
