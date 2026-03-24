"""TotalVI embedding backend.

Requires RNA + protein modalities (Pitfall 3).
"""

from __future__ import annotations

import numpy as np

from sc_tools.assembly.embed._base import register_embedding_backend

__all__ = ["TotalviBackend"]


class TotalviBackend:
    """TotalVI joint embedding backend for RNA + protein data."""

    @staticmethod
    def run(mdata, *, n_factors: int = 15, max_epochs: int = 400, **kwargs) -> tuple[np.ndarray, dict]:
        """Run TotalVI joint embedding.

        Parameters
        ----------
        mdata
            MuData object with RNA and protein modalities.
        n_factors
            Number of latent dimensions.
        max_epochs
            Maximum training epochs.
        **kwargs
            Additional arguments passed to TotalVI.

        Returns
        -------
        embedding
            Array of shape (n_cells, n_factors).
        metadata
            Dict with method parameters.

        Raises
        ------
        SCToolsDataError
            If MuData does not contain both RNA and protein modalities.
        """
        from sc_tools.errors import SCToolsDataError

        mod_keys = set(mdata.mod.keys())
        if "rna" not in mod_keys or "protein" not in mod_keys:
            raise SCToolsDataError(
                "TotalVI requires RNA and protein modalities. "
                f"Found: {sorted(mod_keys)}. "
                "For other combinations, use MOFA+ (default).",
                suggestion="Use --method mofa for non RNA+protein data",
            )

        import scvi  # lazy import per CLI-08

        scvi.model.TOTALVI.setup_mudata(mdata, modalities={"rna": "rna", "protein": "protein"})
        model = scvi.model.TOTALVI(mdata, n_latent=n_factors)
        model.train(max_epochs=max_epochs)
        embedding = model.get_latent_representation()

        obsm_key = "X_totalvi"
        mdata.obsm[obsm_key] = embedding

        metadata = {
            "method": "totalvi",
            "n_factors": n_factors,
            "max_epochs": max_epochs,
        }
        return embedding, metadata


register_embedding_backend("totalvi", TotalviBackend)
