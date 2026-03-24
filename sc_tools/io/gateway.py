"""IO Gateway: tiered data loading with pre-load memory guard.

Tiers
-----
- **T1 (metadata):** h5py-only read -- returns a dict, never loads arrays.
- **T2 (summary):** Backed AnnData -- obs/var accessible, X stays on disk.
- **T3 (full):** Full in-memory AnnData -- guarded by memory check.
"""

from __future__ import annotations

import logging
from enum import Enum
from pathlib import Path

from sc_tools.errors import SCToolsRuntimeError
from sc_tools.io.estimate import estimate_from_h5
from sc_tools.io.metadata import read_h5ad_metadata

logger = logging.getLogger(__name__)


class DataTier(str, Enum):
    """Data loading tier."""

    T1_METADATA = "metadata"
    T2_SUMMARY = "summary"
    T3_FULL = "full"


class IOGateway:
    """Dispatch reads by tier with a memory safety guard on full loads."""

    def read(
        self,
        path: str | Path,
        tier: DataTier,
        *,
        force: bool = False,
    ):
        """Read an h5ad file at the requested tier.

        Parameters
        ----------
        path
            Path to an ``.h5ad`` file.
        tier
            Data loading tier (T1_METADATA, T2_SUMMARY, T3_FULL).
        force
            If ``True``, bypass the memory guard for T3 loads.

        Returns
        -------
        dict | AnnData
            T1 returns a metadata dict. T2/T3 return AnnData.

        Raises
        ------
        SCToolsRuntimeError
            If T3 load would exceed 80% of available RAM (and ``force=False``).
        """
        path = str(path)

        if tier == DataTier.T1_METADATA:
            return read_h5ad_metadata(path)

        if tier == DataTier.T2_SUMMARY:
            import anndata

            return anndata.read_h5ad(path, backed="r")

        # T3_FULL
        if not force:
            self._check_memory_guard(path)

        import anndata

        return anndata.read_h5ad(path)

    def _check_memory_guard(self, path: str | Path) -> None:
        """Raise if estimated peak memory exceeds 80% of available RAM."""
        import psutil

        est = estimate_from_h5(path)
        estimated_peak_mb = est["estimated_peak_mb"]

        available_mb = psutil.virtual_memory().available / (1024**2)

        if estimated_peak_mb > available_mb * 0.8:
            raise SCToolsRuntimeError(
                f"Estimated memory ({estimated_peak_mb:.0f}MB) exceeds 80% of "
                f"available RAM ({available_mb:.0f}MB). "
                "Use --force to override or reduce dataset size."
            )
