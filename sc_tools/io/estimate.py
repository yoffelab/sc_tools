"""Pre-load memory and runtime estimation from h5ad file metadata.

Estimates are computed from h5py metadata (no full data load) so they
can be used for dry-run checks and the ``sct estimate`` command.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from sc_tools.io.metadata import read_h5ad_metadata

# Method-specific multiplier dict (used by Plan 02 estimate command).
METHOD_MULTIPLIERS: dict[str, float] = {
    "preprocess": 2.5,
    "de": 2.0,
    "qc": 1.5,
    "default": 2.0,
}


def estimate_from_h5(path: str | Path) -> dict:
    """Estimate memory requirements for loading an h5ad file.

    Parameters
    ----------
    path
        Path to an ``.h5ad`` file.

    Returns
    -------
    dict
        Keys: ``n_obs``, ``n_vars``, ``x_dtype``, ``x_sparse``,
        ``x_mb``, ``obs_mb``, ``var_mb``, ``base_mb``, ``estimated_peak_mb``.
    """
    meta = read_h5ad_metadata(path)

    n_obs = meta["n_obs"]
    n_vars = meta["n_vars"]
    x_dtype = meta["x_dtype"]
    x_sparse = meta["x_sparse"]

    dtype_size = np.dtype(x_dtype).itemsize

    # --- X memory estimate ---
    if x_sparse:
        nnz = meta.get("nnz", 0)
        # data + indices (int32) + indptr (int32)
        x_bytes = nnz * dtype_size + nnz * 4 + (n_obs + 1) * 4
    else:
        x_bytes = n_obs * n_vars * dtype_size

    # --- obs / var estimates (conservative per-cell/gene overhead) ---
    obs_bytes = n_obs * 200  # ~200 bytes per cell (conservative)
    var_bytes = n_vars * 100  # ~100 bytes per gene

    x_mb = x_bytes / (1024**2)
    obs_mb = obs_bytes / (1024**2)
    var_mb = var_bytes / (1024**2)
    base_mb = x_mb + obs_mb + var_mb

    # Overhead multiplier for processing intermediates
    estimated_peak_mb = base_mb * METHOD_MULTIPLIERS["default"]

    return {
        "n_obs": n_obs,
        "n_vars": n_vars,
        "x_dtype": x_dtype,
        "x_sparse": x_sparse,
        "x_mb": x_mb,
        "obs_mb": obs_mb,
        "var_mb": var_mb,
        "base_mb": base_mb,
        "estimated_peak_mb": estimated_peak_mb,
    }
