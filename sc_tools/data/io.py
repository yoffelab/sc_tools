"""
Data I/O and caching utilities.
"""

import hashlib
import logging
import os
import pickle
from datetime import datetime
from pathlib import Path

logger = logging.getLogger(__name__)


def get_cache_key(
    sc_data_file: str,
    celltype_key: str,
    sc_batch_key: str,
    n_genes_max: int,
    skip_hvg: bool | None = None,
) -> str:
    """
    Generate a cache key based on input parameters and file modification time.

    Parameters
    ----------
    sc_data_file : str
        Path to single-cell data file
    celltype_key : str
        Cell type annotation key
    sc_batch_key : str
        Single-cell batch key
    n_genes_max : int
        Maximum number of genes
    skip_hvg : bool, optional
        Whether to skip HVG computation (for cache key generation)

    Returns
    -------
    str
        MD5 hash of the cache key
    """
    # Include file modification time to detect changes
    try:
        mtime = os.path.getmtime(sc_data_file)
    except OSError:
        mtime = 0

    # Create a hash from parameters
    skip_hvg_str = str(skip_hvg) if skip_hvg is not None else "None"
    key_string = (
        f"{sc_data_file}:{mtime}:{celltype_key}:{sc_batch_key}:{n_genes_max}:{skip_hvg_str}"
    )
    key_hash = hashlib.md5(key_string.encode()).hexdigest()
    return key_hash


def load_cached_signatures(
    cache_path: Path, logger_instance: logging.Logger | None = None
) -> list[str] | None:
    """
    Load cached signature genes from file.

    Parameters
    ----------
    cache_path : Path
        Path to cache file
    logger_instance : Logger, optional
        Custom logger instance. If None, uses module logger.

    Returns
    -------
    list of str or None
        List of signature genes if cache exists and is valid, None otherwise
    """
    log = logger_instance if logger_instance is not None else logger
    try:
        if cache_path.exists():
            with open(cache_path, "rb") as f:
                cache_data = pickle.load(f)
                if isinstance(cache_data, dict) and "signature_genes" in cache_data:
                    log.info(
                        f"   ✅ Loaded {len(cache_data['signature_genes'])} signature genes from cache"
                    )
                    return cache_data["signature_genes"]
                elif isinstance(cache_data, list):
                    # Backward compatibility: old format was just a list
                    log.info(
                        f"   ✅ Loaded {len(cache_data)} signature genes from cache (old format)"
                    )
                    return cache_data
    except Exception as e:
        log.warning(f"   ⚠️  Failed to load cache: {e}")
    return None


def save_cached_signatures(
    signature_genes: list[str],
    cache_path: Path,
    metadata: dict,
    logger_instance: logging.Logger | None = None,
) -> None:
    """
    Save signature genes to cache file with metadata.

    Parameters
    ----------
    signature_genes : list of str
        List of signature gene names
    cache_path : Path
        Path to save cache file
    metadata : dict
        Metadata dictionary to store with cache
    logger_instance : Logger, optional
        Custom logger instance. If None, uses module logger.
    """
    log = logger_instance if logger_instance is not None else logger
    try:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_data = {
            "signature_genes": signature_genes,
            "metadata": metadata,
            "timestamp": datetime.now().isoformat(),
        }
        with open(cache_path, "wb") as f:
            pickle.dump(cache_data, f)
        log.info(f"   💾 Saved {len(signature_genes)} signature genes to cache: {cache_path}")
    except Exception as e:
        log.warning(f"   ⚠️  Failed to save cache: {e}")


def write_h5ad(
    adata,
    basename: str,
    output_dir: str | Path,
    versioned: bool = True,
    dt: datetime | None = None,
) -> Path:
    """
    Save AnnData to h5ad. By default uses versioned filenames for traceability.

    When versioned=True (default): saves to output_dir/YYDDMM.hh.mm.basename.h5ad.
    When versioned=False: saves to output_dir/basename.h5ad (adds .h5ad if missing).

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object to save.
    basename : str
        Base name (no extension), e.g. "adata_genescores".
    output_dir : str or Path
        Directory to write the file; created if missing.
    dt : datetime, optional
        Timestamp for version prefix; if None, uses now().

    Returns
    -------
    Path
        Path to the saved file.
    """
    from ..utils.save import versioned_path

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if versioned:
        out_path = versioned_path(output_dir, basename, "h5ad", subdir=None, dt=dt)
    else:
        name = basename if basename.endswith(".h5ad") else f"{basename}.h5ad"
        out_path = output_dir / name
    adata.write_h5ad(out_path)
    logger.info(f"Saved h5ad: {out_path}")
    return out_path
