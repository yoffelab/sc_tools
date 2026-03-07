"""Runtime environment introspection."""

from __future__ import annotations

import importlib.metadata
import sys


def _pkg_version(name: str) -> str | None:
    """Return installed version of *name*, or None if not installed."""
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def _gpu_devices() -> list[dict] | None:
    """Return list of CUDA devices with name and memory, or None if unavailable."""
    try:
        import torch

        if not torch.cuda.is_available():
            return None
        devices = []
        for i in range(torch.cuda.device_count()):
            props = torch.cuda.get_device_properties(i)
            devices.append(
                {
                    "index": i,
                    "name": props.name,
                    "memory_gb": round(props.total_memory / 1e9, 1),
                }
            )
        return devices or None
    except ImportError:
        return None


def get_environment_info() -> dict:
    """Collect runtime environment information.

    Returns a snapshot of installed package versions and available
    compute backends (GPU, rapids, dask).  All fields are None when
    a package is not installed so callers can adapt without raising.

    Returns
    -------
    dict
        Keys: python, packages (dict of name->version|None),
        gpu (list of device dicts | None), rapids (bool), dask (bool).
    """
    packages = {
        "scanpy": _pkg_version("scanpy"),
        "anndata": _pkg_version("anndata"),
        "scvi-tools": _pkg_version("scvi-tools"),
        "harmonypy": _pkg_version("harmonypy"),
        "rapids-singlecell": _pkg_version("rapids-singlecell"),
        "dask": _pkg_version("dask"),
        "torch": _pkg_version("torch"),
        "numpy": _pkg_version("numpy"),
        "scipy": _pkg_version("scipy"),
        "pandas": _pkg_version("pandas"),
        "squidpy": _pkg_version("squidpy"),
        "spatialdata": _pkg_version("spatialdata"),
    }

    rapids_available = False
    try:
        import rapids_singlecell  # noqa: F401

        rapids_available = True
    except ImportError:
        pass

    dask_available = _pkg_version("dask") is not None

    return {
        "python": sys.version,
        "packages": packages,
        "gpu": _gpu_devices(),
        "rapids": rapids_available,
        "dask": dask_available,
    }
