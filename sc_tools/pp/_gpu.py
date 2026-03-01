"""GPU backend dispatcher for preprocessing.

Auto-detects rapids-singlecell (GPU) and falls back to scanpy (CPU).
All pp functions call get_backend() so the API is identical regardless of hardware.
"""

from __future__ import annotations

import logging
from types import ModuleType

logger = logging.getLogger(__name__)

_backend_cache: tuple[ModuleType, str] | None = None


def get_backend() -> tuple[ModuleType, str]:
    """Return (module, name) for the fastest available backend.

    Returns
    -------
    tuple[ModuleType, str]
        ``(rapids_singlecell, "rapids")`` if available, else ``(scanpy, "scanpy")``.
    """
    global _backend_cache
    if _backend_cache is not None:
        return _backend_cache

    try:
        import rapids_singlecell as rsc

        _backend_cache = (rsc, "rapids")
        logger.info("Using rapids-singlecell (GPU) backend for preprocessing")
    except ImportError:
        import scanpy as sc

        _backend_cache = (sc, "scanpy")
        logger.debug("Using scanpy (CPU) backend for preprocessing")

    return _backend_cache
