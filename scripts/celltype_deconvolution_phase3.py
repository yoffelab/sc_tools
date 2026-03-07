"""
Phase III: Cell Type Deconvolution for Spatial Spots (Modular)

This script performs cell type deconvolution using multiple methods:
- Tangram: Optimal transport-based mapping
- Cell2location: Negative binomial regression model
- DestVI: Bayesian model for cell type proportions

Memory-optimized approach:
- Processes one library_id at a time to avoid memory issues
- Clears variables after each batch
- Uses minimal data copies
- Handles errors gracefully with fallback to alternative methods

Strategy:
- Tries methods in order: DestVI -> Cell2location -> Tangram
- Falls back to next method if one fails
- Extracts cell type proportions from method-specific outputs
- Stores proportions in adata.obsm['cell_type_proportions'] as a matrix
- Also stores individual cell type proportions as obs columns
"""

import gc
import hashlib
import logging
import os
import pickle
import sys
import warnings
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from tqdm import tqdm

# Memory profiling
try:
    import tracemalloc

    import psutil

    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False
    try:
        import tracemalloc
    except ImportError:
        tracemalloc = None

# Suppress warnings for cleaner output
warnings.filterwarnings("ignore", category=UserWarning)

# Setup logging for memory tracking
log_dir = "output/deconvolution/logs"
os.makedirs(log_dir, exist_ok=True)
log_file = os.path.join(log_dir, f"memory_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)

sc.settings.set_figure_params(dpi=300, dpi_save=400)

# ============================================================================
# Configuration
# ============================================================================
CONFIG = {
    "sc_data_file": "results/seurat_object.h5ad",
    "spatial_data_file": (
        "results/adata.scored.h5ad"
        if __import__("os").path.exists("results/adata.scored.h5ad")
        else "results/adata.normalized.scored.p35.h5ad"
    ),
    "celltype_key": "cell.type",
    "spatial_batch_key": "library_id",
    "sc_batch_key": "Batch",
    "output_dir": "output/deconvolution",
    "n_genes_max": 2000,
    "methods": ["destvi", "cell2location", "tangram"],  # Order to try
    "tangram_epochs": 500,
    "cell2location_max_epochs": 25000,
    "destvi_max_epochs": 25000,
    "cache_signature_genes": True,  # Cache signature genes to avoid recomputation
    "signature_cache_dir": "output/deconvolution/cache",  # Directory for cache files
    "force_recompute_signatures": False,  # Force recomputation even if cache exists
}

# ============================================================================
# Memory Profiling Functions
# ============================================================================


def get_memory_usage() -> dict[str, float]:
    """Get current memory usage in MB."""
    memory_info = {}

    if PSUTIL_AVAILABLE:
        process = psutil.Process(os.getpid())
        mem_info = process.memory_info()
        memory_info["rss_mb"] = mem_info.rss / 1024 / 1024  # Resident Set Size
        memory_info["vms_mb"] = mem_info.vms / 1024 / 1024  # Virtual Memory Size
        memory_info["percent"] = process.memory_percent()

        # System memory
        sys_mem = psutil.virtual_memory()
        memory_info["system_available_mb"] = sys_mem.available / 1024 / 1024
        memory_info["system_percent"] = sys_mem.percent
    else:
        memory_info["rss_mb"] = 0.0
        memory_info["vms_mb"] = 0.0
        memory_info["percent"] = 0.0
        memory_info["system_available_mb"] = 0.0
        memory_info["system_percent"] = 0.0

    if tracemalloc is not None and tracemalloc.is_tracing():
        current, peak = tracemalloc.get_traced_memory()
        memory_info["tracemalloc_current_mb"] = current / 1024 / 1024
        memory_info["tracemalloc_peak_mb"] = peak / 1024 / 1024

    return memory_info


def log_memory(step_name: str, adata: ad.AnnData | None = None):
    """Log memory usage at a specific step."""
    mem = get_memory_usage()

    msg = f"[MEMORY] {step_name}:"
    msg += f" RSS={mem['rss_mb']:.1f}MB"
    if mem["percent"] > 0:
        msg += f" ({mem['percent']:.1f}% of process)"
    if mem["system_available_mb"] > 0:
        msg += f" | System: {mem['system_available_mb']:.1f}MB available ({mem['system_percent']:.1f}% used)"
    if "tracemalloc_peak_mb" in mem:
        msg += f" | Peak traced: {mem['tracemalloc_peak_mb']:.1f}MB"

    if adata is not None:
        # Estimate AnnData memory
        if hasattr(adata.X, "data"):
            x_mem = adata.X.data.nbytes / 1024 / 1024
        else:
            x_mem = adata.X.nbytes / 1024 / 1024
        msg += f" | AnnData X: {x_mem:.1f}MB ({adata.shape[0]} spots x {adata.shape[1]} genes)"

    logger.info(msg)
    return mem


def check_memory_threshold(threshold_mb: float = 8000, threshold_percent: float = 85.0) -> bool:
    """Check if memory usage exceeds thresholds."""
    mem = get_memory_usage()

    if mem["rss_mb"] > threshold_mb:
        logger.warning(
            f"[MEMORY WARNING] RSS ({mem['rss_mb']:.1f}MB) exceeds threshold ({threshold_mb}MB)"
        )
        return False

    if mem["system_percent"] > threshold_percent:
        logger.warning(
            f"[MEMORY WARNING] System memory ({mem['system_percent']:.1f}%) exceeds threshold ({threshold_percent}%)"
        )
        return False

    return True


def aggressive_cleanup():
    """Aggressive memory cleanup."""
    gc.collect()
    gc.collect()  # Call twice to handle circular references
    if PSUTIL_AVAILABLE:
        # Force Python to release memory
        import ctypes

        libc = ctypes.CDLL("libc.so.6")
        libc.malloc_trim(0)
    logger.debug("[MEMORY] Aggressive cleanup performed")


def estimate_adata_memory(adata: ad.AnnData) -> float:
    """Estimate memory usage of AnnData object in MB."""
    total = 0

    # X matrix
    if hasattr(adata.X, "data"):
        total += adata.X.data.nbytes
        if hasattr(adata.X, "indices"):
            total += adata.X.indices.nbytes
        if hasattr(adata.X, "indptr"):
            total += adata.X.indptr.nbytes
    else:
        total += adata.X.nbytes

    # obs and var
    total += adata.obs.memory_usage(deep=True).sum()
    total += adata.var.memory_usage(deep=True).sum()

    # obsm
    for key, value in adata.obsm.items():
        if hasattr(value, "nbytes"):
            total += value.nbytes
        elif hasattr(value, "memory_usage"):
            total += value.memory_usage(deep=True).sum()

    return total / 1024 / 1024  # Convert to MB


# ============================================================================
# Data Loading and Preprocessing Functions
# ============================================================================


def load_single_cell_data(
    sc_data_file: str, celltype_key: str, sc_batch_key: str, qc_labels: list[str]
) -> tuple[ad.AnnData, list[str]]:
    """Load and preprocess single-cell reference data."""
    logger.info("=" * 70)
    logger.info("STEP 1: Loading single-cell reference data")
    logger.info("=" * 70)
    log_memory("Before loading single-cell data")

    logger.info(f"Loading from: {sc_data_file}")
    sc_adata = sc.read_h5ad(sc_data_file)
    sc_adata.var_names_make_unique()
    log_memory("After loading single-cell data", sc_adata)
    logger.info(f"Single-cell data: {sc_adata.shape} (cells x genes)")
    logger.info(f"Estimated memory: {estimate_adata_memory(sc_adata):.1f}MB")

    # Preprocess
    log_memory("Before preprocessing")
    if np.max(sc_adata.X) > 100:  # Likely raw counts
        logger.info("Normalizing and log-transforming scRNA-seq data...")
        sc.pp.normalize_total(sc_adata, target_sum=1e4)
        sc.pp.log1p(sc_adata)
        log_memory("After normalization", sc_adata)
    else:
        logger.info("Single-cell data already appears log-normalized")

    # Skip QC filtering for efficiency (as requested)
    logger.info("Skipping QC filtering step for efficiency")
    logger.info(f"Total cells: {sc_adata.n_obs}")

    # Get cell type names (exclude QC labels from cell type list)
    if celltype_key in sc_adata.obs.columns:
        unique_cts = sorted(sc_adata.obs[celltype_key].unique())
        # Filter out QC labels from cell type list, but keep all cells in data
        unique_cts = [ct for ct in unique_cts if ct not in qc_labels]
        logger.info(f"Cell types (excluding QC labels): {len(unique_cts)}")
        logger.info(f"Total unique labels in data: {len(sc_adata.obs[celltype_key].unique())}")
    else:
        raise ValueError(f"celltype_key '{celltype_key}' not found in single-cell data")

    log_memory("Final single-cell data state", sc_adata)
    return sc_adata, unique_cts


def get_cache_key(sc_data_file: str, celltype_key: str, sc_batch_key: str, n_genes_max: int) -> str:
    """Generate a cache key based on input parameters and file modification time."""
    # Include file modification time to detect changes
    try:
        mtime = os.path.getmtime(sc_data_file)
    except OSError:
        mtime = 0

    # Create a hash from parameters
    key_string = f"{sc_data_file}:{mtime}:{celltype_key}:{sc_batch_key}:{n_genes_max}"
    key_hash = hashlib.md5(key_string.encode()).hexdigest()
    return key_hash


def load_cached_signatures(cache_path: Path) -> list[str] | None:
    """Load cached signature genes from file."""
    try:
        if cache_path.exists():
            with open(cache_path, "rb") as f:
                cache_data = pickle.load(f)
                if isinstance(cache_data, dict) and "signature_genes" in cache_data:
                    logger.info(
                        f"   ✅ Loaded {len(cache_data['signature_genes'])} signature genes from cache"
                    )
                    return cache_data["signature_genes"]
                elif isinstance(cache_data, list):
                    # Backward compatibility: old format was just a list
                    logger.info(
                        f"   ✅ Loaded {len(cache_data)} signature genes from cache (old format)"
                    )
                    return cache_data
    except Exception as e:
        logger.warning(f"   ⚠️  Failed to load cache: {e}")
    return None


def save_cached_signatures(signature_genes: list[str], cache_path: Path, metadata: dict) -> None:
    """Save signature genes to cache file with metadata."""
    try:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_data = {
            "signature_genes": signature_genes,
            "metadata": metadata,
            "timestamp": datetime.now().isoformat(),
        }
        with open(cache_path, "wb") as f:
            pickle.dump(cache_data, f)
        logger.info(f"   💾 Saved {len(signature_genes)} signature genes to cache: {cache_path}")
    except Exception as e:
        logger.warning(f"   ⚠️  Failed to save cache: {e}")


def select_signature_genes(
    sc_adata: ad.AnnData,
    celltype_key: str,
    sc_batch_key: str,
    n_genes_max: int,
    cache_dir: str | None = None,
    force_recompute: bool = False,
    sc_data_file: str | None = None,
) -> list[str]:
    """Select signature genes for deconvolution."""
    logger.info("\n2. Selecting signature genes for deconvolution...")

    # Check cache if enabled
    if cache_dir and not force_recompute and sc_data_file:
        cache_key = get_cache_key(sc_data_file, celltype_key, sc_batch_key, n_genes_max)
        cache_path = Path(cache_dir) / f"signature_genes_{cache_key}.pkl"

        cached_genes = load_cached_signatures(cache_path)
        if cached_genes is not None:
            return cached_genes
        logger.info("   Cache miss or invalid, computing signature genes...")
    elif force_recompute:
        logger.info("   Force recompute enabled, computing signature genes...")

    # Step 1: Highly variable genes
    logger.info("   Computing highly variable genes...")
    sc.pp.highly_variable_genes(
        sc_adata, flavor="seurat_v3", n_top_genes=3000, subset=False, batch_key=sc_batch_key
    )
    hvg_genes = set(sc_adata.var[sc_adata.var.highly_variable].index)
    logger.info(f"   Found {len(hvg_genes)} HVGs")

    # Step 2: Top marker genes per cell type
    logger.info("   Computing marker genes per cell type...")
    sc.tl.rank_genes_groups(sc_adata, groupby=celltype_key, method="wilcoxon", use_raw=False)
    marker_genes = []
    for group in sc_adata.obs[celltype_key].unique():
        df = sc.get.rank_genes_groups_df(sc_adata, group=group)
        top_genes = df.head(100).names.tolist()
        marker_genes.extend(top_genes)
    marker_genes = set(marker_genes)
    logger.info(f"   Found {len(marker_genes)} marker genes")

    # Step 3: Combine
    candidate_genes = hvg_genes.union(marker_genes)
    logger.info(f"   Candidate genes: {len(candidate_genes)}")

    candidate_genes_list = list(candidate_genes)

    # Save to cache if enabled
    if cache_dir and sc_data_file:
        cache_key = get_cache_key(sc_data_file, celltype_key, sc_batch_key, n_genes_max)
        cache_path = Path(cache_dir) / f"signature_genes_{cache_key}.pkl"
        metadata = {
            "sc_data_file": sc_data_file,
            "celltype_key": celltype_key,
            "sc_batch_key": sc_batch_key,
            "n_genes_max": n_genes_max,
            "n_genes": len(candidate_genes_list),
        }
        save_cached_signatures(candidate_genes_list, cache_path, metadata)

    return candidate_genes_list


def get_library_ids(spatial_data_file: str, spatial_batch_key: str) -> list[str]:
    """Get library IDs from spatial data without loading full dataset."""
    print("\n3. Identifying libraries to process...")
    spatial_metadata = pd.read_h5ad(spatial_data_file, backed="r").obs
    if spatial_batch_key not in spatial_metadata.columns:
        if "sample_id" in spatial_metadata.columns:
            spatial_batch_key = "sample_id"
            print("   Using 'sample_id' instead of 'library_id'")
        else:
            raise ValueError("Neither 'library_id' nor 'sample_id' found in spatial data")

    library_ids = sorted(spatial_metadata[spatial_batch_key].unique())
    print(f"   Found {len(library_ids)} libraries to process: {library_ids}")
    return library_ids, spatial_batch_key


# ============================================================================
# Deconvolution Methods
# ============================================================================


def deconvolve_tangram(
    sc_adata: ad.AnnData,
    spatial_adata_lib: ad.AnnData,
    shared_genes: list[str],
    celltype_key: str,
    output_path: str,
    num_epochs: int,
) -> np.ndarray | None:
    """Run Tangram deconvolution for a single library."""
    try:
        import tangram as tg

        log_memory("Tangram: Start")

        # Check if already computed
        if os.path.exists(output_path):
            logger.info("Loading existing Tangram result...")
            ad_map = sc.read_h5ad(output_path)
            log_memory("Tangram: After loading existing result")
        else:
            logger.info("Running Tangram mapping...")
            log_memory("Tangram: Before creating copies")

            # Create minimal copies with only needed genes
            sc_copy = sc_adata[:, shared_genes].copy()
            spatial_lib_sig = spatial_adata_lib[:, shared_genes].copy()
            log_memory("Tangram: After creating gene subsets", sc_copy)
            log_memory("Tangram: Spatial subset", spatial_lib_sig)

            # Check memory before Tangram
            if not check_memory_threshold(threshold_mb=12000, threshold_percent=90.0):
                logger.warning("[MEMORY] Memory too high for Tangram, skipping...")
                sc_copy = None
                spatial_lib_sig = None
                aggressive_cleanup()
                return None

            logger.info("Tangram: Preparing data...")
            tg.pp_adatas(sc_copy, spatial_lib_sig, genes=shared_genes)
            log_memory("Tangram: After pp_adatas")

            logger.info("Tangram: Running map_cells_to_space (this may take time)...")
            ad_map = tg.map_cells_to_space(
                sc_copy,
                spatial_lib_sig,
                mode="clusters",
                cluster_label=celltype_key,
                num_epochs=num_epochs,
            )
            log_memory("Tangram: After map_cells_to_space")

            # Clear copies before saving
            sc_copy = None
            spatial_lib_sig = None
            aggressive_cleanup()

            ad_map.write(output_path)
            logger.info(f"Saved: {output_path}")
            log_memory("Tangram: After saving")

        # Project annotations
        tg.project_cell_annotations(ad_map, spatial_adata_lib, annotation=celltype_key)

        # Extract proportions
        unique_cts = sorted(sc_adata.obs[celltype_key].unique())
        proportions = None

        # Try obsm first
        for key in spatial_adata_lib.obsm.keys():
            if "proportion" in key.lower() or celltype_key.lower() in key.lower():
                prop_data = spatial_adata_lib.obsm[key]
                if prop_data.shape[0] == spatial_adata_lib.n_obs:
                    proportions = prop_data
                    break

        # Try obs columns
        if proportions is None:
            ct_cols = [col for col in spatial_adata_lib.obs.columns if col in unique_cts]
            if len(ct_cols) >= len(unique_cts) * 0.8:
                proportions_list = []
                for ct in unique_cts:
                    if ct in spatial_adata_lib.obs.columns:
                        proportions_list.append(spatial_adata_lib.obs[ct].values)
                    else:
                        proportions_list.append(np.zeros(spatial_adata_lib.n_obs))
                proportions = np.column_stack(proportions_list)

        return proportions

    except Exception as e:
        logger.error(f"ERROR with Tangram: {e}")
        import traceback

        logger.error(traceback.format_exc())
        log_memory("Tangram: After error")
        aggressive_cleanup()
        return None


def deconvolve_cell2location(
    sc_adata: ad.AnnData,
    spatial_adata_lib: ad.AnnData,
    shared_genes: list[str],
    celltype_key: str,
    output_path: str,
    max_epochs: int,
) -> np.ndarray | None:
    """Run Cell2location deconvolution for a single library."""
    try:
        import cell2location as c2l
        from cell2location.models import Cell2location, RegressionModel

        logger.info("Running Cell2location...")
        log_memory("Cell2location: Start")

        # Check memory before starting
        if not check_memory_threshold(threshold_mb=12000, threshold_percent=90.0):
            logger.warning("[MEMORY] Memory too high for Cell2location, skipping...")
            return None

        # Prepare data - Cell2location needs raw counts
        logger.info("Cell2location: Creating gene subsets...")
        sc_copy = sc_adata[:, shared_genes].copy()
        spatial_lib_sig = spatial_adata_lib[:, shared_genes].copy()
        log_memory("Cell2location: After creating subsets", sc_copy)

        # Ensure raw counts (minimal copies)
        if sc_copy.raw is None:
            sc_copy.raw = sc_copy.copy()
        if spatial_lib_sig.raw is None:
            spatial_lib_sig.raw = spatial_lib_sig.copy()
        log_memory("Cell2location: After setting raw")

        # Train regression model on single-cell data
        logger.info("Cell2location: Training regression model on single-cell data...")
        c2l.models.RegressionModel.setup_anndata(sc_copy, labels_key=celltype_key, batch_key=None)
        mod = RegressionModel(sc_copy)
        log_memory("Cell2location: After creating regression model")

        # Reduce epochs/samples if memory is tight
        actual_epochs = max_epochs
        num_samples = 1000
        batch_size = 2500

        if not check_memory_threshold(threshold_mb=10000, threshold_percent=85.0):
            actual_epochs = min(max_epochs, 15000)
            num_samples = 500
            batch_size = 1500
            logger.info(
                f"Cell2location: Reducing parameters (epochs={actual_epochs}, samples={num_samples}) due to memory"
            )

        mod.train(max_epochs=actual_epochs, use_gpu=False)
        log_memory("Cell2location: After training regression model")

        sc_copy = mod.export_posterior(
            sc_copy, sample_kwargs={"num_samples": num_samples, "batch_size": batch_size}
        )
        log_memory("Cell2location: After exporting posterior")

        # Train Cell2location on spatial data
        logger.info("Cell2location: Training on spatial data...")
        c2l.models.Cell2location.setup_anndata(spatial_lib_sig, batch_key=None)
        mod_spatial = Cell2location(spatial_lib_sig, sc_copy)
        log_memory("Cell2location: After creating spatial model")

        mod_spatial.train(max_epochs=actual_epochs, use_gpu=False)
        log_memory("Cell2location: After training spatial model")

        # Get proportions
        logger.info("Cell2location: Extracting proportions...")
        spatial_lib_sig = mod_spatial.export_posterior(
            spatial_lib_sig, sample_kwargs={"num_samples": num_samples, "batch_size": batch_size}
        )
        log_memory("Cell2location: After exporting spatial posterior")

        # Extract proportions from obsm
        # Cell2location stores cell abundance in multiple quantiles
        # Use mean or q05 (5th percentile) for conservative estimates
        if "q05_cell_abundance_w_sf" in spatial_lib_sig.obsm:
            proportions = spatial_lib_sig.obsm["q05_cell_abundance_w_sf"].values.copy()
        elif "means_cell_abundance_w_sf" in spatial_lib_sig.obsm:
            proportions = spatial_lib_sig.obsm["means_cell_abundance_w_sf"].values.copy()
        else:
            logger.warning("Could not find Cell2location proportions in obsm")
            logger.warning(f"Available obsm keys: {list(spatial_lib_sig.obsm.keys())}")
            mod = None
            mod_spatial = None
            sc_copy = None
            spatial_lib_sig = None
            aggressive_cleanup()
            return None

        # Normalize to proportions (cell abundance -> proportions)
        row_sums = proportions.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1  # Avoid division by zero
        proportions = proportions / row_sums

        # Cleanup
        mod = None
        mod_spatial = None
        sc_copy = None
        spatial_lib_sig = None
        aggressive_cleanup()
        log_memory("Cell2location: After cleanup")

        return proportions

    except ImportError:
        logger.warning("Cell2location not installed, skipping...")
        return None
    except Exception as e:
        logger.error(f"ERROR with Cell2location: {e}")
        import traceback

        logger.error(traceback.format_exc())
        log_memory("Cell2location: After error")
        aggressive_cleanup()
        return None


def deconvolve_destvi(
    sc_adata: ad.AnnData,
    spatial_adata_lib: ad.AnnData,
    shared_genes: list[str],
    celltype_key: str,
    output_path: str,
    max_epochs: int,
) -> np.ndarray | None:
    """Run DestVI deconvolution for a single library."""
    try:
        import scvi
        from scvi.external import DestVI

        logger.info("Running DestVI...")
        log_memory("DestVI: Start")

        # Check memory before starting
        if not check_memory_threshold(threshold_mb=12000, threshold_percent=90.0):
            logger.warning("[MEMORY] Memory too high for DestVI, skipping...")
            return None

        # Prepare data - use minimal copies
        logger.info("DestVI: Creating gene subsets...")
        sc_copy = sc_adata[:, shared_genes].copy()
        spatial_lib_sig = spatial_adata_lib[:, shared_genes].copy()
        log_memory("DestVI: After creating subsets", sc_copy)

        # Ensure raw counts (but avoid unnecessary copies)
        if sc_copy.raw is None:
            sc_copy.raw = sc_copy.copy()
        if spatial_lib_sig.raw is None:
            spatial_lib_sig.raw = spatial_lib_sig.copy()
        log_memory("DestVI: After setting raw")

        # Setup and train on single-cell data
        logger.info("DestVI: Training on single-cell data...")
        DestVI.setup_anndata(sc_copy, labels_key=celltype_key)
        vae = DestVI(sc_copy)
        log_memory("DestVI: After creating VAE model")

        # Reduce epochs if memory is tight
        actual_epochs = max_epochs
        if not check_memory_threshold(threshold_mb=10000, threshold_percent=85.0):
            actual_epochs = min(max_epochs, 10000)
            logger.info(f"DestVI: Reducing epochs to {actual_epochs} due to memory constraints")

        vae.train(max_epochs=actual_epochs)
        log_memory("DestVI: After training VAE")

        # Clear single-cell copy before spatial training
        sc_copy = None
        aggressive_cleanup()

        # Train on spatial data
        logger.info("DestVI: Training on spatial data...")
        DestVI.setup_anndata(spatial_lib_sig, labels_key=None)
        vae_st = DestVI.load_query_data(spatial_lib_sig, vae)
        log_memory("DestVI: After loading query data")

        vae_st.train(max_epochs=actual_epochs)
        log_memory("DestVI: After training spatial model")

        # Get proportions
        logger.info("DestVI: Extracting proportions...")
        spatial_lib_sig.obsm["proportions"] = vae_st.get_proportions(spatial_lib_sig)
        proportions = spatial_lib_sig.obsm["proportions"].values.copy()

        # Cleanup
        vae = None
        vae_st = None
        spatial_lib_sig = None
        aggressive_cleanup()
        log_memory("DestVI: After cleanup")

        return proportions

    except ImportError:
        logger.warning("DestVI not installed, skipping...")
        return None
    except Exception as e:
        logger.error(f"ERROR with DestVI: {e}")
        import traceback

        logger.error(traceback.format_exc())
        log_memory("DestVI: After error")
        aggressive_cleanup()
        return None


# ============================================================================
# Main Processing Function
# ============================================================================


def process_library(
    lib_id: str,
    sc_adata: ad.AnnData,
    spatial_data_file: str,
    spatial_batch_key: str,
    candidate_genes: list[str],
    celltype_key: str,
    unique_cts: list[str],
    methods: list[str],
    config: dict,
) -> tuple[np.ndarray | None, list[str]]:
    """Process a single library with fallback methods."""
    logger.info(f"\n{'=' * 60}")
    logger.info(f"Processing Library: {lib_id}")
    logger.info(f"{'=' * 60}")
    log_memory(f"Start processing {lib_id}")

    try:
        # Check memory before loading
        if not check_memory_threshold():
            logger.warning(
                f"[MEMORY] Memory usage high before loading {lib_id}, performing cleanup..."
            )
            aggressive_cleanup()

        # Load library data
        logger.info(f"Loading spatial data for {lib_id}...")
        spatial_adata_full = sc.read_h5ad(spatial_data_file, backed="r")
        log_memory(f"After opening backed file for {lib_id}")

        lib_mask = spatial_adata_full.obs[spatial_batch_key] == lib_id
        lib_indices = spatial_adata_full.obs[lib_mask].index

        if len(lib_indices) == 0:
            logger.warning(f"No spots found for library {lib_id}, skipping")
            spatial_adata_full = None
            aggressive_cleanup()
            return None, []

        logger.info(f"Loading {len(lib_indices)} spots for {lib_id}...")
        spatial_adata_lib = spatial_adata_full[lib_mask].copy()
        spatial_adata_full = None
        aggressive_cleanup()
        log_memory(f"After loading {lib_id} data", spatial_adata_lib)

        logger.info(f"Library {lib_id}: {spatial_adata_lib.shape} (spots x genes)")
        logger.info(f"Estimated memory: {estimate_adata_memory(spatial_adata_lib):.1f}MB")

        # Use raw counts if available
        if spatial_adata_lib.raw is not None:
            logger.info("Using raw counts from spatial data")
            log_memory("Before loading raw data")
            raw_adata = spatial_adata_lib.raw.to_adata()
            log_memory("After loading raw data", raw_adata)
            spatial_adata_lib.X = raw_adata[:, spatial_adata_lib.var_names].X
            raw_adata = None
            aggressive_cleanup()
            log_memory("After replacing X with raw counts", spatial_adata_lib)

        # Find shared genes
        shared_genes = list(
            set(candidate_genes).intersection(sc_adata.var_names, spatial_adata_lib.var_names)
        )[: config["n_genes_max"]]
        logger.info(f"Shared signature genes: {len(shared_genes)}")

        if len(shared_genes) < 100:
            logger.warning(f"Too few shared genes ({len(shared_genes)}), skipping")
            spatial_adata_lib = None
            aggressive_cleanup()
            return None, []

        # Check memory before processing
        if not check_memory_threshold():
            logger.warning(
                f"[MEMORY] Memory usage high before processing {lib_id}, reducing gene set..."
            )
            # Reduce number of genes if memory is tight
            if len(shared_genes) > 1000:
                shared_genes = shared_genes[:1000]
                logger.info(f"Reduced to {len(shared_genes)} genes due to memory constraints")
            aggressive_cleanup()

        # Try each method in order
        proportions = None
        successful_methods = []

        for method in methods:
            logger.info(f"\nTrying method: {method.upper()}")
            log_memory(f"Before {method} for {lib_id}", spatial_adata_lib)

            # Check memory before each method
            if not check_memory_threshold(threshold_mb=10000, threshold_percent=90.0):
                logger.warning(
                    f"[MEMORY] Memory usage too high for {method}, skipping to next method..."
                )
                aggressive_cleanup()
                continue

            output_path = os.path.join(config["output_dir"], f"{method}_pred_{lib_id}.h5ad")

            if method == "tangram":
                proportions = deconvolve_tangram(
                    sc_adata,
                    spatial_adata_lib,
                    shared_genes,
                    celltype_key,
                    output_path,
                    config["tangram_epochs"],
                )
            elif method == "cell2location":
                proportions = deconvolve_cell2location(
                    sc_adata,
                    spatial_adata_lib,
                    shared_genes,
                    celltype_key,
                    output_path,
                    config["cell2location_max_epochs"],
                )
            elif method == "destvi":
                proportions = deconvolve_destvi(
                    sc_adata,
                    spatial_adata_lib,
                    shared_genes,
                    celltype_key,
                    output_path,
                    config["destvi_max_epochs"],
                )
            else:
                print(f"   Unknown method: {method}, skipping")
                continue

            if proportions is not None:
                # Validate proportions shape
                if proportions.shape[0] == spatial_adata_lib.n_obs:
                    # Check if we need to align cell types
                    if proportions.shape[1] != len(unique_cts):
                        print(
                            f"   Warning: Proportions shape {proportions.shape} doesn't match expected ({spatial_adata_lib.n_obs}, {len(unique_cts)})"
                        )
                        # Try to match cell types if possible
                        # For now, accept if number of spots matches
                        if proportions.shape[1] > 0:
                            # Pad or truncate to match expected number of cell types
                            if proportions.shape[1] < len(unique_cts):
                                # Pad with zeros
                                padding = np.zeros(
                                    (proportions.shape[0], len(unique_cts) - proportions.shape[1])
                                )
                                proportions = np.hstack([proportions, padding])
                            else:
                                # Truncate
                                proportions = proportions[:, : len(unique_cts)]
                            # Renormalize
                            row_sums = proportions.sum(axis=1, keepdims=True)
                            row_sums[row_sums == 0] = 1
                            proportions = proportions / row_sums

                    successful_methods.append(method)
                    logger.info(f"✓ Successfully processed {lib_id} with {method.upper()}")
                    log_memory(f"After successful {method} for {lib_id}")
                    break
                else:
                    logger.warning(
                        f"Proportions shape mismatch ({proportions.shape[0]} spots vs {spatial_adata_lib.n_obs}), trying next method..."
                    )
                    proportions = None
                    aggressive_cleanup()

            # Cleanup after each method attempt
            log_memory(f"After {method} attempt for {lib_id}")
            aggressive_cleanup()

        # Final cleanup
        spatial_adata_lib = None
        aggressive_cleanup()
        log_memory(f"After processing {lib_id}")

        return proportions, successful_methods

    except Exception as e:
        logger.error(f"ERROR processing library {lib_id}: {e}")
        import traceback

        logger.error(traceback.format_exc())
        log_memory(f"After error in {lib_id}")
        aggressive_cleanup()
        return None, []


def integrate_proportions(
    all_proportions: list[np.ndarray],
    all_spot_indices: list[str],
    spatial_data_file: str,
    unique_cts: list[str],
    celltype_key: str,
) -> None:
    """Integrate proportions from all libraries into main spatial data."""
    print("\n5. Integrating cell type proportions into spatial data...")

    if len(all_proportions) == 0:
        print("   ⚠️  WARNING: No cell type proportions extracted!")
        return

    # Concatenate proportions
    proportions_matrix = np.vstack(all_proportions)
    print(f"   Total proportions matrix: {proportions_matrix.shape} (spots x cell types)")

    # Load full spatial data
    print("   Loading full spatial data for integration...")
    spatial_adata = sc.read_h5ad(spatial_data_file)

    # Reorder proportions to match spatial_adata order
    if len(all_spot_indices) == proportions_matrix.shape[0]:
        index_map = {idx: i for i, idx in enumerate(all_spot_indices)}
        ordered_indices = [index_map.get(idx, -1) for idx in spatial_adata.obs_names]

        valid_mask = np.array(ordered_indices) >= 0
        ordered_proportions = np.zeros((spatial_adata.n_obs, proportions_matrix.shape[1]))
        ordered_proportions[valid_mask] = proportions_matrix[np.array(ordered_indices)[valid_mask]]

        # Store in obsm
        spatial_adata.obsm["cell_type_proportions"] = ordered_proportions

        # Store as individual obs columns
        for i, ct in enumerate(unique_cts):
            if i < ordered_proportions.shape[1]:
                col_name = f"prop_{ct.replace(' ', '_').replace('/', '_')}"
                spatial_adata.obs[col_name] = ordered_proportions[:, i]

        print(f"   Stored proportions in obsm and {len(unique_cts)} obs columns")

        # Save integrated data
        output_file = "results/adata.deconvolution.h5ad"
        spatial_adata.write(output_file)
        print(f"\n   ✓ Saved integrated data: {output_file}")

        # Summary statistics
        print("\n   Cell type proportion summary:")
        mean_props = ordered_proportions.mean(axis=0)
        for ct, prop in zip(unique_cts[: len(mean_props)], mean_props):
            print(f"      {ct}: {prop:.4f}")


# ============================================================================
# Main Execution
# ============================================================================


def main():
    """Main execution function."""
    logger.info("=" * 70)
    logger.info("PHASE III: CELL TYPE DECONVOLUTION (Modular, Memory-Optimized)")
    logger.info("=" * 70)

    # Start memory tracking
    if tracemalloc is not None:
        tracemalloc.start()
        logger.info("Started tracemalloc for memory tracking")

    log_memory("Script start")

    # Create output directories
    os.makedirs(CONFIG["output_dir"], exist_ok=True)
    os.makedirs("figures/deconvolution", exist_ok=True)

    # QC labels
    qc_labels = ["QC_Filtered", "Doublets", "Low quality", "Unknown III (SM)"]

    # Load single-cell data
    sc_adata, unique_cts = load_single_cell_data(
        CONFIG["sc_data_file"], CONFIG["celltype_key"], CONFIG["sc_batch_key"], qc_labels
    )

    # Select signature genes
    candidate_genes = select_signature_genes(
        sc_adata,
        CONFIG["celltype_key"],
        CONFIG["sc_batch_key"],
        CONFIG["n_genes_max"],
        cache_dir=CONFIG.get("signature_cache_dir")
        if CONFIG.get("cache_signature_genes", False)
        else None,
        force_recompute=CONFIG.get("force_recompute_signatures", False),
        sc_data_file=CONFIG["sc_data_file"],
    )

    # Get library IDs
    library_ids, spatial_batch_key = get_library_ids(
        CONFIG["spatial_data_file"], CONFIG["spatial_batch_key"]
    )

    # Process each library
    print("\n4. Processing libraries one at a time (memory-efficient)...")
    print(f"   Methods to try (in order): {CONFIG['methods']}")

    all_proportions = []
    all_spot_indices = []
    processed_libraries = []
    method_usage = dict.fromkeys(CONFIG["methods"], 0)

    for lib_id in tqdm(library_ids, desc="Processing libraries"):
        proportions, successful_methods = process_library(
            lib_id,
            sc_adata,
            CONFIG["spatial_data_file"],
            spatial_batch_key,
            candidate_genes,
            CONFIG["celltype_key"],
            unique_cts,
            CONFIG["methods"],
            CONFIG,
        )

        if proportions is not None:
            all_proportions.append(proportions)
            # Get spot indices for this library
            spatial_metadata = pd.read_h5ad(CONFIG["spatial_data_file"], backed="r").obs
            lib_mask = spatial_metadata[spatial_batch_key] == lib_id
            all_spot_indices.extend(spatial_metadata[lib_mask].index.tolist())
            processed_libraries.append(lib_id)

            for method in successful_methods:
                method_usage[method] += 1

    # Integrate proportions
    integrate_proportions(
        all_proportions,
        all_spot_indices,
        CONFIG["spatial_data_file"],
        unique_cts,
        CONFIG["celltype_key"],
    )

    # Summary
    log_memory("Script end")
    logger.info("\n" + "=" * 70)
    logger.info("✅ CELL TYPE DECONVOLUTION COMPLETE")
    logger.info("=" * 70)
    logger.info(f"Processed libraries: {len(processed_libraries)}/{len(library_ids)}")
    logger.info("Method usage:")
    for method, count in method_usage.items():
        logger.info(f"   {method.upper()}: {count} libraries")
    if len(processed_libraries) > 0:
        logger.info("Output file: results/adata.deconvolution.h5ad")
        logger.info("Cell type proportions stored in: adata.obsm['cell_type_proportions']")
    logger.info(f"Memory log saved to: {log_file}")

    if tracemalloc is not None:
        current, peak = tracemalloc.get_traced_memory()
        logger.info(f"Tracemalloc peak memory: {peak / 1024 / 1024:.1f}MB")
        tracemalloc.stop()


if __name__ == "__main__":
    main()
