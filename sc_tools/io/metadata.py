"""Read h5ad metadata via h5py without loading array data into memory.

The fallback chain for n_obs/n_vars mirrors the pattern used in
``sc_tools.mcp.tools_server.inspect_checkpoint``.
"""

from __future__ import annotations

from pathlib import Path


def read_h5ad_metadata(path: str | Path) -> dict:
    """Extract metadata from an h5ad file using h5py (no full data load).

    Parameters
    ----------
    path
        Path to an ``.h5ad`` file.

    Returns
    -------
    dict
        Keys: ``n_obs``, ``n_vars``, ``x_dtype``, ``x_sparse``,
        ``obs_columns``, ``obsm_keys``, ``layer_keys``, ``uns_keys``.
        If sparse, also includes ``nnz``.
    """
    import h5py

    path = str(path)

    with h5py.File(path, "r") as f:
        # --- n_obs via fallback chain ---
        n_obs = f["obs"].attrs.get("_index_length", None)
        if n_obs is None and "_index" in f["obs"]:
            n_obs = len(f["obs"]["_index"])
        elif n_obs is None:
            first_key = next(iter(f["obs"].keys()), None)
            n_obs = len(f["obs"][first_key]) if first_key else 0
        n_obs = int(n_obs)

        # --- n_vars via fallback chain ---
        n_vars = f["var"].attrs.get("_index_length", None)
        if n_vars is None and "_index" in f["var"]:
            n_vars = len(f["var"]["_index"])
        elif n_vars is None:
            first_key = next(iter(f["var"].keys()), None)
            n_vars = len(f["var"][first_key]) if first_key else 0
        n_vars = int(n_vars)

        # --- X dtype and sparsity ---
        x_sparse = isinstance(f["X"], h5py.Group)
        if x_sparse:
            x_dtype = str(f["X"]["data"].dtype)
            nnz = int(f["X"]["data"].shape[0])
        else:
            x_dtype = str(f["X"].dtype)
            nnz = None

        # --- Structured metadata ---
        obs_columns = [k for k in f["obs"].keys() if not k.startswith("__")]
        obsm_keys = list(f["obsm"].keys()) if "obsm" in f else []
        layer_keys = list(f["layers"].keys()) if "layers" in f else []
        uns_keys = list(f["uns"].keys()) if "uns" in f else []

    result: dict = {
        "n_obs": n_obs,
        "n_vars": n_vars,
        "x_dtype": x_dtype,
        "x_sparse": x_sparse,
        "obs_columns": obs_columns,
        "obsm_keys": obsm_keys,
        "layer_keys": layer_keys,
        "uns_keys": uns_keys,
    }
    if nnz is not None:
        result["nnz"] = nnz

    return result
