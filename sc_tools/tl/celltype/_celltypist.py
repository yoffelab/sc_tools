"""
CellTypist backend (soft dependency: celltypist>=1.3).

Install: ``pip install 'sc-tools[celltyping]'``
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import anndata as ad

from ._base import register_celltype_backend


class CellTypistBackend:
    """CellTypist automated cell typing; requires celltypist soft dep."""

    @staticmethod
    def run(
        adata: ad.AnnData,
        *,
        cluster_key: str = "leiden",
        store_proba: bool = False,
        model: str = "Immune_All_Low.pkl",
        majority_voting: bool = True,
        p_thres: float = 0.5,
        **kwargs,
    ) -> tuple[pd.Series, pd.Series, pd.DataFrame | None, dict]:
        """
        Run CellTypist annotation.

        Parameters
        ----------
        adata
            AnnData; expression should be log1p-normalised counts in ``X``.
        cluster_key
            Obs key for cluster labels (used for majority voting).
        store_proba
            Store per-celltype probability DataFrame in obsm.
        model
            CellTypist model name or path.
        majority_voting
            Apply majority voting per cluster after initial annotation.
        p_thres
            Probability threshold for assignment.

        Returns
        -------
        labels, scores, proba_or_None, metadata
        """
        try:
            import celltypist
        except (ImportError, TypeError) as err:
            raise ImportError(
                "celltypist is required for the 'celltypist' backend. "
                "pip install 'sc-tools[celltyping]' or pip install celltypist"
            ) from err

        predictions = celltypist.annotate(
            adata,
            model=model,
            majority_voting=majority_voting,
            over_clustering=cluster_key if majority_voting else None,
            p_thres=p_thres,
        )

        result_adata = predictions.to_adata()
        label_col = "majority_voting" if majority_voting else "predicted_labels"
        labels = result_adata.obs[label_col].reindex(adata.obs_names)

        # Confidence scores from conf_score column
        if "conf_score" in result_adata.obs.columns:
            scores = result_adata.obs["conf_score"].reindex(adata.obs_names).astype(np.float64)
        else:
            scores = pd.Series(np.ones(adata.n_obs, dtype=np.float64), index=adata.obs_names)

        # Probability matrix
        proba: pd.DataFrame | None = None
        if (
            store_proba
            and hasattr(result_adata, "obsm")
            and "celltypist_prob_matrix" in result_adata.obsm
        ):
            proba_arr = result_adata.obsm["celltypist_prob_matrix"]
            proba = pd.DataFrame(proba_arr, index=adata.obs_names)

        meta = {
            "model": model,
            "majority_voting": majority_voting,
            "p_thres": p_thres,
        }
        return labels, scores, proba, meta


register_celltype_backend("celltypist", CellTypistBackend)
