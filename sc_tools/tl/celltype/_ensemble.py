"""
Ensemble voting backend (no soft dependencies).

Reads previously stored ``adata.obs[f'celltype_auto_{m}']`` columns for each
method in ``methods`` and combines them into a consensus label.

Strategies
----------
majority_vote
    Each cell is assigned the most frequent label across methods.
    Ties broken by alphabetical order (deterministic).
weighted_vote
    Each method is weighted by its confidence score
    (``adata.obs[f'celltype_auto_{m}_score']``). The label with the
    highest sum of weights wins.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import anndata as ad

from ._base import register_celltype_backend


class EnsembleBackend:
    """Ensemble voting over previously run annotation methods."""

    @staticmethod
    def run(
        adata: ad.AnnData,
        *,
        cluster_key: str = "leiden",
        store_proba: bool = False,
        methods: list[str],
        strategy: str = "majority_vote",
        **kwargs,
    ) -> tuple[pd.Series, pd.Series, pd.DataFrame | None, dict]:
        """
        Combine multiple celltype_auto_* columns by voting.

        Parameters
        ----------
        adata
            AnnData with prior annotation columns in obs.
        cluster_key
            Unused; present for API compatibility.
        store_proba
            Not applicable; always returns None.
        methods
            List of method names (short names, not full obs key). Each must
            have a corresponding ``celltype_auto_{m}`` column in obs.
        strategy
            ``'majority_vote'`` or ``'weighted_vote'``.

        Returns
        -------
        labels, scores, None, metadata
        """
        if not methods:
            raise ValueError("methods list cannot be empty for ensemble voting.")

        label_cols: list[str] = []
        score_cols: list[str] = []
        for m in methods:
            lk = f"celltype_auto_{m}"
            sk = f"celltype_auto_{m}_score"
            if lk not in adata.obs.columns:
                raise KeyError(
                    f"Column '{lk}' not found in adata.obs. "
                    f"Run annotate_celltypes(adata, method='{m}') first."
                )
            label_cols.append(lk)
            score_cols.append(sk)

        label_matrix = adata.obs[label_cols].astype(str)  # (n_obs, n_methods)
        n_obs = adata.n_obs

        if strategy == "majority_vote":
            labels_arr = np.full(n_obs, "Unknown", dtype=object)
            scores_arr = np.zeros(n_obs, dtype=np.float64)
            for i in range(n_obs):
                row = label_matrix.iloc[i].values
                unique, counts = np.unique(row, return_counts=True)
                best_idx = np.argmax(counts)
                best_count = counts[best_idx]
                # If tie, pick alphabetically first
                ties = unique[counts == best_count]
                labels_arr[i] = sorted(ties)[0]
                scores_arr[i] = best_count / len(methods)

        elif strategy == "weighted_vote":
            labels_arr = np.full(n_obs, "Unknown", dtype=object)
            scores_arr = np.zeros(n_obs, dtype=np.float64)
            # Collect weights; use uniform 1.0 if score column missing
            weight_list = []
            for sk in score_cols:
                if sk in adata.obs.columns:
                    weight_list.append(adata.obs[sk].values.astype(np.float64))
                else:
                    weight_list.append(np.ones(n_obs, dtype=np.float64))
            weight_matrix = np.stack(weight_list, axis=1)  # (n_obs, n_methods)

            for i in range(n_obs):
                row_labels = label_matrix.iloc[i].values
                row_weights = weight_matrix[i]
                # Sum weights per label
                agg: dict[str, float] = {}
                for lbl, w in zip(row_labels, row_weights, strict=True):
                    agg[lbl] = agg.get(lbl, 0.0) + float(w)
                best_lbl = max(agg, key=lambda x: (agg[x], x))  # stable tie-break
                labels_arr[i] = best_lbl
                total_w = sum(row_weights)
                scores_arr[i] = agg[best_lbl] / total_w if total_w > 0 else 0.0

        else:
            raise ValueError(
                f"Unknown ensemble strategy '{strategy}'. Use 'majority_vote' or 'weighted_vote'."
            )

        labels = pd.Series(labels_arr, index=adata.obs_names, dtype=object)
        scores = pd.Series(scores_arr, index=adata.obs_names, dtype=np.float64)
        meta = {"methods": list(methods), "strategy": strategy}
        return labels, scores, None, meta


register_celltype_backend("ensemble", EnsembleBackend)
