 """
Memory profiling and management utilities.

Provides functions for tracking memory usage, performing cleanup,
and estimating memory requirements for AnnData objects.
"""

import os
import sys
import gc
import logging
from typing import Dict, Optional
import anndata as ad

logger = logging.getLogger(__name__)

# Try to import optional dependencies
try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False

try:
    import tracemalloc
    TRACEMALLOC_AVAILABLE = True
except ImportError:
    tracemalloc = None
    TRACEMALLOC_AVAILABLE = False


def get_memory_usage() -> Dict[str, float]:
    """
    Get current memory usage in MB.
    
    Returns
    -------
    dict
        Dictionary with memory usage metrics:
        - 'rss_mb': Resident Set Size in MB
        - 'vms_mb': Virtual Memory Size in MB (if psutil available)
        - 'percent': Memory usage as percentage of process (if psutil available)
        - 'system_available_mb': Available system memory in MB (if psutil available)
        - 'system_percent': System memory usage percentage (if psutil available)
        - 'tracemalloc_current_mb': Current traced memory (if tracemalloc active)
        - 'tracemalloc_peak_mb': Peak traced memory (if tracemalloc active)
    """
    memory_info = {}
    
    if PSUTIL_AVAILABLE:
        process = psutil.Process(os.getpid())
        mem_info = process.memory_info()
        memory_info['rss_mb'] = mem_info.rss / 1024 / 1024  # Resident Set Size
        memory_info['vms_mb'] = mem_info.vms / 1024 / 1024  # Virtual Memory Size
        memory_info['percent'] = process.memory_percent()
        
        # System memory
        sys_mem = psutil.virtual_memory()
        memory_info['system_available_mb'] = sys_mem.available / 1024 / 1024
        memory_info['system_percent'] = sys_mem.percent
    else:
        memory_info['rss_mb'] = 0.0
        memory_info['vms_mb'] = 0.0
        memory_info['percent'] = 0.0
        memory_info['system_available_mb'] = 0.0
        memory_info['system_percent'] = 0.0
    
    if TRACEMALLOC_AVAILABLE and tracemalloc is not None and tracemalloc.is_tracing():
        current, peak = tracemalloc.get_traced_memory()
        memory_info['tracemalloc_current_mb'] = current / 1024 / 1024
        memory_info['tracemalloc_peak_mb'] = peak / 1024 / 1024
    
    return memory_info


def log_memory(step_name: str, adata: Optional[ad.AnnData] = None, 
               logger_instance: Optional[logging.Logger] = None) -> Dict[str, float]:
    """
    Log memory usage at a specific step.
    
    Parameters
    ----------
    step_name : str
        Name of the step for logging
    adata : AnnData, optional
        Optional AnnData object to estimate its memory usage
    logger_instance : Logger, optional
        Custom logger instance. If None, uses module logger.
    
    Returns
    -------
    dict
        Memory usage dictionary (from get_memory_usage)
    """
    log = logger_instance if logger_instance is not None else logger
    mem = get_memory_usage()
    
    msg = f"[MEMORY] {step_name}:"
    msg += f" RSS={mem['rss_mb']:.1f}MB"
    if mem['percent'] > 0:
        msg += f" ({mem['percent']:.1f}% of process)"
    if mem['system_available_mb'] > 0:
        msg += f" | System: {mem['system_available_mb']:.1f}MB available ({mem['system_percent']:.1f}% used)"
    if 'tracemalloc_peak_mb' in mem:
        msg += f" | Peak traced: {mem['tracemalloc_peak_mb']:.1f}MB"
    
    if adata is not None:
        # Estimate AnnData memory
        x_mem = estimate_adata_memory(adata)
        msg += f" | AnnData X: {x_mem:.1f}MB ({adata.shape[0]} spots x {adata.shape[1]} genes)"
    
    log.info(msg)
    return mem


def aggressive_cleanup():
    """
    Aggressively clean up memory.
    
    Performs garbage collection and attempts to release memory
    back to the operating system (platform-dependent).
    """
    gc.collect()
    gc.collect()  # Call twice to handle circular references
    
    if PSUTIL_AVAILABLE:
        # Force Python to release memory (Linux/macOS)
        try:
            import ctypes
            if sys.platform == "darwin":
                libc = ctypes.CDLL("libc.dylib")
            else:
                libc = ctypes.CDLL("libc.so.6")
            libc.malloc_trim(0)
        except (OSError, AttributeError):
            pass  # Not available on this platform


def estimate_adata_memory(adata: ad.AnnData) -> float:
    """
    Estimate memory usage of AnnData object in MB.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object to estimate
    
    Returns
    -------
    float
        Estimated memory usage in MB
    """
    total = 0
    
    # X matrix
    if hasattr(adata.X, 'data'):
        total += adata.X.data.nbytes
        if hasattr(adata.X, 'indices'):
            total += adata.X.indices.nbytes
        if hasattr(adata.X, 'indptr'):
            total += adata.X.indptr.nbytes
    else:
        total += adata.X.nbytes
    
    # obs and var
    total += adata.obs.memory_usage(deep=True).sum()
    total += adata.var.memory_usage(deep=True).sum()
    
    # obsm
    for key, value in adata.obsm.items():
        if hasattr(value, 'nbytes'):
            total += value.nbytes
        elif hasattr(value, 'memory_usage'):
            total += value.memory_usage(deep=True).sum()
    
    return total / 1024 / 1024  # Convert to MB


def check_memory_threshold(threshold_mb: float = 8000, 
                          threshold_percent: float = 85.0,
                          logger_instance: Optional[logging.Logger] = None) -> bool:
    """
    Check if memory usage exceeds thresholds.
    
    Parameters
    ----------
    threshold_mb : float
        Maximum RSS memory in MB
    threshold_percent : float
        Maximum system memory usage percentage
    logger_instance : Logger, optional
        Custom logger instance. If None, uses module logger.
    
    Returns
    -------
    bool
        True if memory is below thresholds, False otherwise
    """
    log = logger_instance if logger_instance is not None else logger
    mem = get_memory_usage()
    
    if mem['rss_mb'] > threshold_mb:
        log.warning(f"[MEMORY WARNING] RSS ({mem['rss_mb']:.1f}MB) exceeds threshold ({threshold_mb}MB)")
        return False
    
    if mem['system_percent'] > threshold_percent:
        log.warning(f"[MEMORY WARNING] System memory ({mem['system_percent']:.1f}%) exceeds threshold ({threshold_percent}%)")
        return False
    
    return True
