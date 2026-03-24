"""
Data I/O and caching utilities.
"""

from __future__ import annotations

import hashlib
import logging
import os
import pickle
from datetime import datetime
from pathlib import Path

import pandas as pd

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


def _coerce_arrow_strings(adata) -> None:
    """Convert Arrow-backed string columns/index to plain Python types.

    Newer pandas (2.x) + pyarrow produce ``ArrowStringArray`` from parquet/CSV
    reads and h5ad round-trips.  The ``h5py`` backend used by
    ``anndata.write_h5ad`` cannot serialize these, raising::

        IORegistryError: No method registered for writing
        <class 'pandas.arrays.ArrowStringArray'> into <class 'h5py._hl.group.Group'>

    This has bitten us repeatedly (run_qc_report, run_preprocessing,
    run_integration_benchmark, robin/robin_utils, robin/bin_categorization,
    robin/run_signature_scoring_v2).  Centralising the fix here means every
    caller of ``write_h5ad`` is protected automatically.

    See Also
    --------
    Known error doc: sc_tools/docs/known-errors/arrow-string-h5ad.md
    """
    for df in [adata.obs, adata.var]:
        for col in df.columns:
            if isinstance(df[col].dtype, pd.CategoricalDtype):
                # Rebuild categories — they may be Arrow-backed underneath
                df[col] = df[col].astype(str).astype("category")
            elif not pd.api.types.is_numeric_dtype(df[col]) and df[col].dtype.name != "bool":
                # Coerce any non-numeric, non-bool column to plain object
                df[col] = df[col].astype(str).astype(object)
        # Index coercion
        if not pd.api.types.is_numeric_dtype(df.index):
            df.index = df.index.astype(str)


def write_h5ad(
    adata,
    basename: str,
    output_dir: str | Path,
    versioned: bool = True,
    dt: datetime | None = None,
) -> Path:
    """
    Save AnnData to h5ad. By default uses versioned filenames for traceability.

    Automatically coerces Arrow-backed string columns to plain Python types
    before writing, preventing ``IORegistryError`` on HPC anndata installs.

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
    _coerce_arrow_strings(adata)
    adata.write_h5ad(out_path)
    logger.info(f"Saved h5ad: {out_path}")
    return out_path
