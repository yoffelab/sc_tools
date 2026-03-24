"""MOFA+ embedding backend.

Uses muon.tl.mofa with use_obs='union' for outer join across modalities (D-01).
Recommended default for multi-omic data (D-07).
"""

from __future__ import annotations

import numpy as np

from sc_tools.assembly.embed._base import register_embedding_backend

__all__ = ["MofaBackend"]


class MofaBackend:
    """MOFA+ joint embedding backend."""

    @staticmethod
    def run(mdata, *, n_factors: int = 15, gpu_mode: bool = False, **kwargs) -> tuple[np.ndarray, dict]:
        """Run MOFA+ joint embedding.

        Parameters
        ----------
        mdata
            MuData object with multi-modal data.
        n_factors
            Number of latent factors.
        gpu_mode
            Whether to use GPU acceleration.
        **kwargs
            Additional arguments passed to mu.tl.mofa.

        Returns
        -------
        embedding
            Array of shape (n_cells, n_factors).
        metadata
            Dict with method parameters.
        """
        import muon as mu  # lazy import per CLI-08

        mu.tl.mofa(
            mdata,
            use_obs="union",
            n_factors=n_factors,
            gpu_mode=gpu_mode,
            **kwargs,
        )

        embedding = mdata.obsm["X_mofa"]
        metadata = {
            "method": "mofa",
            "n_factors": n_factors,
            "gpu_mode": gpu_mode,
        }
        return embedding, metadata


register_embedding_backend("mofa", MofaBackend)
