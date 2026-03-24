"""MultiVI embedding backend.

Requires RNA + ATAC modalities (Pitfall 3).
"""

from __future__ import annotations

import numpy as np

from sc_tools.assembly.embed._base import register_embedding_backend

__all__ = ["MultiviBackend"]


class MultiviBackend:
    """MultiVI joint embedding backend for RNA + ATAC data."""

    @staticmethod
    def run(mdata, *, n_factors: int = 15, max_epochs: int = 400, **kwargs) -> tuple[np.ndarray, dict]:
        """Run MultiVI joint embedding.

        Parameters
        ----------
        mdata
            MuData object with RNA and ATAC modalities.
        n_factors
            Number of latent dimensions.
        max_epochs
            Maximum training epochs.
        **kwargs
            Additional arguments passed to MultiVI.

        Returns
        -------
        embedding
            Array of shape (n_cells, n_factors).
        metadata
            Dict with method parameters.

        Raises
        ------
        SCToolsDataError
            If MuData does not contain both RNA and ATAC modalities.
        """
        from sc_tools.errors import SCToolsDataError

        mod_keys = set(mdata.mod.keys())
        if "rna" not in mod_keys or "atac" not in mod_keys:
            raise SCToolsDataError(
                "MultiVI requires RNA and ATAC modalities. "
                f"Found: {sorted(mod_keys)}. "
                "For other combinations, use MOFA+ (default).",
                suggestion="Use --method mofa for non RNA+ATAC data",
            )

        import scvi  # lazy import per CLI-08

        scvi.model.MULTIVI.setup_mudata(mdata, modalities={"rna": "rna", "atac": "atac"})
        model = scvi.model.MULTIVI(mdata, n_latent=n_factors)
        model.train(max_epochs=max_epochs)
        embedding = model.get_latent_representation()

        obsm_key = "X_multivi"
        mdata.obsm[obsm_key] = embedding

        metadata = {
            "method": "multivi",
            "n_factors": n_factors,
            "max_epochs": max_epochs,
        }
        return embedding, metadata


register_embedding_backend("multivi", MultiviBackend)
