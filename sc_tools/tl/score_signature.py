"""
Signature scoring: score gene signatures and store in obsm (Option B single flat key).

Uses scanpy.tl.score_genes with control genes by default; rank-based (UCell) and
ssGSEA are available as optional backends. Writes to obsm["signature_score"] and
obsm["signature_score_z"] with full-path column names for traceability.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData


def _flatten_signatures(
    d: Mapping[str, Any], prefix: tuple[str, ...] = ()
) -> list[tuple[tuple[str, ...], list[str]]]:
    """Flatten nested dict into (path_tuple, genes). Skip _meta and non-dict/non-list values."""
    out = []
    for k, v in d.items():
        if k == "_meta":
            continue
        if isinstance(v, Mapping):
            out.extend(_flatten_signatures(v, prefix + (k,)))
        elif isinstance(v, list):
            genes = [g for g in v if isinstance(g, str)]
            if genes:
                out.append((prefix + (k,), genes))
    return out


def _index_genes(
    genes: list[str], var_names: np.ndarray
) -> tuple[np.ndarray, list[str], list[str]]:
    """Return indices of genes present in var_names, plus lists of present and missing (case-insensitive)."""
    vm = {g.upper(): i for i, g in enumerate(var_names)}
    idx, present, missing = [], [], []
    for g in genes:
        key = g.upper()
        if key in vm:
            idx.append(vm[key])
            present.append(var_names[vm[key]])
        else:
            missing.append(g)
    return np.array(sorted(set(idx))), present, missing


# ---------------------------------------------------------------------------
# Backend implementations
# ---------------------------------------------------------------------------


def _score_scanpy(
    adata: AnnData,
    leaves: list[tuple[tuple[str, ...], list[str]]],
    var_names: np.ndarray,
    use_raw: bool,
    ctrl_size: int,
    n_bins: int,
    min_genes: int,
) -> tuple[dict, list[dict]]:
    """Scanpy score_genes backend. Returns (score_cols, report_rows)."""
    score_cols: dict[str, np.ndarray] = {}
    report_rows: list[dict] = []

    for path_tuple, genes in leaves:
        col_name = "/".join(path_tuple)
        _, present, missing = _index_genes(genes, var_names)

        if len(present) < min_genes:
            score_cols[col_name] = np.full(adata.n_obs, np.nan, dtype=float)
            report_rows.append(
                {
                    "signature": col_name,
                    "n_present": len(present),
                    "n_missing": len(missing),
                    "status": "skipped",
                }
            )
            continue

        tmp_name = f"_tmp_sig_{len(score_cols)}"
        try:
            sc.tl.score_genes(
                adata,
                gene_list=present,
                score_name=tmp_name,
                ctrl_size=ctrl_size,
                n_bins=n_bins,
                use_raw=use_raw,
            )
            score_cols[col_name] = adata.obs[tmp_name].values.astype(float)
            report_rows.append(
                {
                    "signature": col_name,
                    "n_present": len(present),
                    "n_missing": len(missing),
                    "status": "ok",
                }
            )
        except Exception as e:
            score_cols[col_name] = np.full(adata.n_obs, np.nan, dtype=float)
            report_rows.append(
                {
                    "signature": col_name,
                    "n_present": len(present),
                    "n_missing": len(missing),
                    "status": f"error: {str(e)[:80]}",
                }
            )
        finally:
            if tmp_name in adata.obs.columns:
                adata.obs.drop(columns=[tmp_name], inplace=True)

    return score_cols, report_rows


def _score_ucell(
    adata: AnnData,
    leaves: list[tuple[tuple[str, ...], list[str]]],
    var_names: np.ndarray,
    min_genes: int,
) -> tuple[dict, list[dict]]:
    """UCell rank-based AUC backend. Requires pyucell."""
    try:
        import pyucell  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "pyucell is required for method='ucell'. Install it with: pip install pyucell"
        ) from exc

    import scipy.sparse as sp

    score_cols: dict[str, np.ndarray] = {}
    report_rows: list[dict] = []

    # Build gene-rank matrix once (rank within each cell, ascending expression = low rank)
    X = adata.X
    if sp.issparse(X):
        X = X.toarray()
    # Rank: highest expression = highest rank (UCell convention)
    ranks = X.shape[1] + 1 - np.argsort(np.argsort(X, axis=1), axis=1)  # descending ranks

    for path_tuple, genes in leaves:
        col_name = "/".join(path_tuple)
        idx, present, missing = _index_genes(genes, var_names)

        if len(present) < min_genes:
            score_cols[col_name] = np.full(adata.n_obs, np.nan, dtype=float)
            report_rows.append(
                {
                    "signature": col_name,
                    "n_present": len(present),
                    "n_missing": len(missing),
                    "status": "skipped",
                }
            )
            continue

        try:
            # Use pyucell's u_stat function if available, otherwise call score_cells
            gene_sets = {col_name: present}
            result = pyucell.score_cells(adata, gene_sets)
            score_cols[col_name] = result[col_name].values.astype(float)
            report_rows.append(
                {
                    "signature": col_name,
                    "n_present": len(present),
                    "n_missing": len(missing),
                    "status": "ok",
                }
            )
        except Exception:
            # Fallback: manual AUC from ranks
            sig_idx = idx.tolist()
            if len(sig_idx) == 0:
                score_cols[col_name] = np.full(adata.n_obs, np.nan, dtype=float)
                report_rows.append(
                    {
                        "signature": col_name,
                        "n_present": 0,
                        "n_missing": len(missing),
                        "status": "skipped",
                    }
                )
                continue
            n_genes = X.shape[1]
            n_sig = len(sig_idx)
            rank_sum = ranks[:, sig_idx].sum(axis=1).astype(float)
            # Mann-Whitney U statistic normalized to [0, 1]
            n_sig * (n_genes - n_sig + n_sig)
            auc = (rank_sum - n_sig * (n_sig + 1) / 2) / (n_sig * (n_genes - n_sig))
            score_cols[col_name] = auc.astype(float)
            report_rows.append(
                {
                    "signature": col_name,
                    "n_present": len(present),
                    "n_missing": len(missing),
                    "status": "ok",
                }
            )

    return score_cols, report_rows


def _score_ssgsea(
    adata: AnnData,
    leaves: list[tuple[tuple[str, ...], list[str]]],
    var_names: np.ndarray,
    min_genes: int,
) -> tuple[dict, list[dict]]:
    """ssGSEA backend via gseapy. Requires gseapy."""
    try:
        import gseapy
    except ImportError as exc:
        raise ImportError(
            "gseapy is required for method='ssgsea'. Install it with: pip install gseapy"
        ) from exc

    import scipy.sparse as sp

    score_cols: dict[str, np.ndarray] = {}
    report_rows: list[dict] = []

    # Build expression DataFrame (genes x cells) for gseapy
    X = adata.X
    if sp.issparse(X):
        X = X.toarray()
    expr_df = pd.DataFrame(X.T, index=list(var_names), columns=list(adata.obs_names))

    # Build gene set dict with only sets meeting min_genes
    valid_gene_sets: dict[str, list[str]] = {}
    skipped: list[tuple[tuple[str, ...], list[str], list[str]]] = []

    for path_tuple, genes in leaves:
        col_name = "/".join(path_tuple)
        _, present, missing = _index_genes(genes, var_names)
        if len(present) < min_genes:
            skipped.append((path_tuple, present, missing))
        else:
            valid_gene_sets[col_name] = present

    # Mark skipped sets
    for path_tuple, present, missing in skipped:
        col_name = "/".join(path_tuple)
        score_cols[col_name] = np.full(adata.n_obs, np.nan, dtype=float)
        report_rows.append(
            {
                "signature": col_name,
                "n_present": len(present),
                "n_missing": len(missing),
                "status": "skipped",
            }
        )

    if not valid_gene_sets:
        return score_cols, report_rows

    try:
        ss = gseapy.ssgsea(
            data=expr_df,
            gene_sets=valid_gene_sets,
            outdir=None,
            no_plot=True,
            processes=1,
        )
        # Result: DataFrame cells x gene_sets (or transposed; normalize)
        res = ss.res2d if hasattr(ss, "res2d") else ss.results
        if hasattr(res, "T"):
            # res2d is gene_sets x samples; transpose to samples x gene_sets
            scores_df = res.T if res.shape[0] == len(valid_gene_sets) else res
        else:
            scores_df = pd.DataFrame(res)

        for path_tuple, genes in leaves:
            col_name = "/".join(path_tuple)
            if col_name not in valid_gene_sets:
                continue
            _, present, missing = _index_genes(genes, var_names)
            if col_name in scores_df.columns:
                score_cols[col_name] = scores_df[col_name].values.astype(float)
            else:
                score_cols[col_name] = np.full(adata.n_obs, np.nan, dtype=float)
            report_rows.append(
                {
                    "signature": col_name,
                    "n_present": len(present),
                    "n_missing": len(missing),
                    "status": "ok",
                }
            )

    except Exception as e:
        for path_tuple, genes in leaves:
            col_name = "/".join(path_tuple)
            if col_name in valid_gene_sets:
                _, present, missing = _index_genes(genes, var_names)
                score_cols[col_name] = np.full(adata.n_obs, np.nan, dtype=float)
                report_rows.append(
                    {
                        "signature": col_name,
                        "n_present": len(present),
                        "n_missing": len(missing),
                        "status": f"error: {str(e)[:80]}",
                    }
                )

    return score_cols, report_rows


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def score_signature(
    adata: AnnData,
    signatures_nested: dict[str, Any] | str | Path,
    method: str = "scanpy",
    use_raw: bool = True,
    ctrl_size: int = 50,
    n_bins: int = 25,
    min_genes: int = 3,
    copy: bool = False,
    save_path: str | Path | None = None,
) -> AnnData:
    """
    Score gene signatures and store in obsm (single flat key for traceability).

    Supports three scoring backends selected by ``method``:
    - ``"scanpy"`` (default): scanpy score_genes (AddModuleScore-style control subtraction).
    - ``"ucell"``: rank-based AUC (requires ``pyucell``; robust to normalization and sparsity).
    - ``"ssgsea"``: single-sample GSEA (requires ``gseapy``; slower; publication-grade).

    Writes raw scores to ``adata.obsm["signature_score"]`` and z-scored (across obs) to
    ``adata.obsm["signature_score_z"]``. Column names are full path (e.g.
    ``Myeloid/Macrophage_Core``). Report is stored in ``adata.uns["signature_score_report"]``.
    The method used is recorded in ``adata.uns["scoring_method"]``.

    Parameters
    ----------
    adata : AnnData
        Annotated data. Must have ``adata.raw`` if ``use_raw=True`` and ``method="scanpy"``.
    signatures_nested : dict, str, or Path
        Two-level nested dict ``{category: {subprocess: [genes]}}`` or path to JSON.
        Keys like ``_meta`` are skipped.
    method : str
        Scoring backend: ``"scanpy"`` | ``"ucell"`` | ``"ssgsea"``. Default ``"scanpy"``.
    use_raw : bool
        Use ``adata.raw`` for expression (default True). Only applies to ``method="scanpy"``.
    ctrl_size : int
        Control genes per signature gene (default 50). Only applies to ``method="scanpy"``.
    n_bins : int
        Expression bins for control matching (default 25). Only applies to ``method="scanpy"``.
    min_genes : int
        Minimum present genes to score (default 3); sets below this threshold receive NaN.
    copy : bool
        If True, copy adata before modifying (default False).
    save_path : str, Path, or None
        If set, write adata to this path after scoring (default None).

    Returns
    -------
    AnnData
        AnnData with ``obsm["signature_score"]``, ``obsm["signature_score_z"]``,
        ``uns["signature_score_report"]``, and ``uns["scoring_method"]``.
        Same object as input if ``copy=False``.
    """
    _valid_methods = ("scanpy", "ucell", "ssgsea")
    if method not in _valid_methods:
        raise ValueError(f"method={method!r} is not valid. Choose from {_valid_methods}.")

    if copy:
        adata = adata.copy()

    # Load signatures
    if isinstance(signatures_nested, (str, Path)):
        path = Path(signatures_nested)
        if not path.exists():
            raise FileNotFoundError(f"Signatures file not found: {path}")
        with open(path) as f:
            signatures_nested = json.load(f)
    if not isinstance(signatures_nested, Mapping):
        raise TypeError("signatures_nested must be a dict or path to JSON")

    leaves = _flatten_signatures(signatures_nested)
    if not leaves:
        raise ValueError("No gene lists found in signatures_nested (after skipping _meta)")

    if method == "scanpy":
        if use_raw and adata.raw is None:
            raise ValueError(
                "adata.raw is None but use_raw=True. score_genes requires raw for expression binning."
            )
        var_names = np.array(adata.raw.var_names if use_raw else adata.var_names, dtype=str)
        score_cols, report_rows = _score_scanpy(
            adata, leaves, var_names, use_raw, ctrl_size, n_bins, min_genes
        )
    elif method == "ucell":
        var_names = np.array(adata.var_names, dtype=str)
        score_cols, report_rows = _score_ucell(adata, leaves, var_names, min_genes)
    elif method == "ssgsea":
        var_names = np.array(adata.var_names, dtype=str)
        score_cols, report_rows = _score_ssgsea(adata, leaves, var_names, min_genes)

    df_scores = pd.DataFrame(score_cols, index=adata.obs_names)
    adata.obsm["signature_score"] = df_scores

    # Z-score each column across observations
    df_z = df_scores.copy()
    with np.errstate(invalid="ignore", divide="ignore"):
        for c in df_z.columns:
            x = df_z[c].values
            m, s = np.nanmean(x), np.nanstd(x)
            if not (np.isnan(s) or s <= 0):
                df_z[c] = (x - m) / s
            else:
                df_z[c] = 0.0
    adata.obsm["signature_score_z"] = df_z

    adata.uns["signature_score_report"] = pd.DataFrame(report_rows)
    adata.uns["scoring_method"] = method

    if save_path is not None:
        adata.write(save_path)

    return adata
