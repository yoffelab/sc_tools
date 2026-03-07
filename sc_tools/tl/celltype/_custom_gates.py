"""
Hierarchical threshold gating backend (zero soft dependencies).

Gate format (dict or YAML file path)::

    {
        "Tcell": {
            "CD3E": (0.5, None),    # (min_expr, max_expr), None = unbounded
            "CD19": (None, 0.1),
            "children": {
                "CD4_Tcell": {"CD4": (0.5, None)},
                "CD8_Tcell": {"CD8A": (0.5, None)},
            }
        }
    }

DFS traversal; each cell is assigned the deepest matching gate.
Cells that do not match any gate receive the label ``"Unassigned"``.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
import scipy.sparse as sp

if TYPE_CHECKING:
    import anndata as ad

from ._base import register_celltype_backend


def _load_gates(gates: dict | str | Path) -> dict:
    """Load gate spec from dict or YAML/JSON file."""
    if isinstance(gates, dict):
        return gates
    path = Path(gates)
    if path.suffix in (".yaml", ".yml"):
        try:
            import yaml
        except ImportError as e:
            raise ImportError("PyYAML required for YAML gate files. pip install pyyaml") from e
        with open(path) as fh:
            return yaml.safe_load(fh)
    import json

    with open(path) as fh:
        return json.load(fh)


def _cell_passes_gate(
    cell_expr: np.ndarray,
    var_upper: dict[str, int],
    gate_spec: dict[str, Any],
) -> bool:
    """Check whether a single cell passes all marker thresholds in gate_spec (ignoring 'children')."""
    for marker, bounds in gate_spec.items():
        if marker == "children":
            continue
        if not isinstance(bounds, (list, tuple)) or len(bounds) != 2:
            continue
        lo, hi = bounds
        idx = var_upper.get(marker.upper())
        if idx is None:
            # Marker not in adata — treat as missing, gate fails
            return False
        val = float(cell_expr[idx])
        if lo is not None and val < lo:
            return False
        if hi is not None and val > hi:
            return False
    return True


def _assign_cell_dfs(
    cell_expr: np.ndarray,
    var_upper: dict[str, int],
    gates: dict[str, Any],
    parent_label: str | None = None,
) -> str:
    """
    Recursively traverse gates DFS and return deepest matching label.

    Returns "Unassigned" if no top-level gate is satisfied.
    """
    for label, spec in gates.items():
        if not _cell_passes_gate(cell_expr, var_upper, spec):
            continue
        # This cell passes the current gate; try children
        children = spec.get("children", {})
        if children:
            child_label = _assign_cell_dfs(cell_expr, var_upper, children, parent_label=label)
            if child_label != "Unassigned":
                return child_label
        # No child matched (or no children) — return current label
        return label
    return "Unassigned"


class CustomGatesBackend:
    """Hierarchical threshold gating; pure numpy, zero soft deps."""

    @staticmethod
    def run(
        adata: ad.AnnData,
        *,
        cluster_key: str = "leiden",
        store_proba: bool = False,
        gates: dict | str | Path,
        normalize: bool = False,
        **kwargs,
    ) -> tuple[pd.Series, pd.Series, pd.DataFrame | None, dict]:
        """
        Apply hierarchical threshold gating.

        Parameters
        ----------
        adata
            AnnData; expression in ``X``.
        cluster_key
            Obs key for cluster labels (used for metadata only; gating is per-cell).
        store_proba
            Not applicable for gating; always returns None.
        gates
            Gate specification dict or path to YAML/JSON file.
        normalize
            If True, normalise each cell to sum to 1 before applying thresholds.

        Returns
        -------
        labels, scores, None, metadata
        """
        gate_spec = _load_gates(gates)

        # Dense expression matrix
        X_src = adata.X
        if sp.issparse(X_src):
            X_dense = X_src.toarray().astype(np.float32)
        else:
            X_dense = np.asarray(X_src, dtype=np.float32)

        if normalize:
            row_sums = X_dense.sum(axis=1, keepdims=True)
            row_sums[row_sums == 0] = 1.0
            X_dense = X_dense / row_sums

        var_upper = {g.upper(): i for i, g in enumerate(adata.var_names)}

        n_obs = adata.n_obs
        labels_arr = np.full(n_obs, "Unassigned", dtype=object)
        scores_arr = np.zeros(n_obs, dtype=np.float64)

        for i in range(n_obs):
            label = _assign_cell_dfs(X_dense[i], var_upper, gate_spec)
            labels_arr[i] = label
            scores_arr[i] = 1.0 if label != "Unassigned" else 0.0

        labels = pd.Series(labels_arr, index=adata.obs_names, dtype=object)
        scores = pd.Series(scores_arr, index=adata.obs_names, dtype=np.float64)

        n_top_gates = len([k for k in gate_spec if not k.startswith("_")])
        meta = {"n_top_gates": n_top_gates, "normalize": normalize}
        return labels, scores, None, meta


register_celltype_backend("custom_gates", CustomGatesBackend)
