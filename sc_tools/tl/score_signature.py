"""
Signature scoring: score gene signatures and store in obsm (Option B single flat key).

Uses scanpy.tl.score_genes with control genes; writes to obsm["signature_score"]
and obsm["signature_score_z"] with full-path column names for traceability.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping, Union

import numpy as np
import pandas as pd
import scanpy as sc

from anndata import AnnData


def _flatten_signatures(d: Mapping[str, Any], prefix: tuple[str, ...] = ()) -> list[tuple[tuple[str, ...], list[str]]]:
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


def _index_genes(genes: list[str], var_names: np.ndarray) -> tuple[np.ndarray, list[str], list[str]]:
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


def score_signature(
    adata: AnnData,
    signatures_nested: Union[dict[str, Any], str, Path],
    use_raw: bool = True,
    ctrl_size: int = 50,
    n_bins: int = 25,
    min_genes: int = 3,
    copy: bool = False,
    save_path: Union[str, Path, None] = None,
) -> AnnData:
    """
    Score gene signatures and store in obsm (single flat key for traceability).

    Uses scanpy.tl.score_genes with control genes. Writes raw scores to
    adata.obsm["signature_score"] and z-scored (across obs) to
    adata.obsm["signature_score_z"]. Column names are full path (e.g.
    Myeloid/Macrophage_Core). Report is stored in adata.uns["signature_score_report"].

    Parameters
    ----------
    adata : AnnData
        Annotated data. Must have adata.raw if use_raw=True.
    signatures_nested : dict, str, or Path
        Two-level nested dict {category: {subprocess: [genes]}} or path to JSON.
        Keys like _meta are skipped.
    use_raw : bool
        Use adata.raw for expression (default True). score_genes needs raw for binning.
    ctrl_size : int
        Control genes per signature gene (default 50).
    n_bins : int
        Expression bins for control matching (default 25).
    min_genes : int
        Minimum present genes to score (default 3); else NaN and status skipped.
    copy : bool
        If True, copy adata before modifying (default False).
    save_path : str, Path, or None
        If set, write adata to this path after scoring (default None).

    Returns
    -------
    AnnData
        AnnData with obsm["signature_score"], obsm["signature_score_z"], and
        uns["signature_score_report"]. Same object as input if copy=False.
    """
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

    if use_raw and adata.raw is None:
        raise ValueError("adata.raw is None but use_raw=True. score_genes requires raw for expression binning.")

    var_names = np.array(adata.raw.var_names if use_raw else adata.var_names, dtype=str)

    score_cols = {}
    report_rows = []

    for path_tuple, genes in leaves:
        col_name = "/".join(path_tuple)
        idx, present, missing = _index_genes(genes, var_names)

        if len(present) < min_genes:
            score_cols[col_name] = np.full(adata.n_obs, np.nan, dtype=float)
            report_rows.append({
                "signature": col_name,
                "n_present": len(present),
                "n_missing": len(missing),
                "status": "skipped",
            })
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
            report_rows.append({
                "signature": col_name,
                "n_present": len(present),
                "n_missing": len(missing),
                "status": "ok",
            })
        except Exception as e:
            score_cols[col_name] = np.full(adata.n_obs, np.nan, dtype=float)
            report_rows.append({
                "signature": col_name,
                "n_present": len(present),
                "n_missing": len(missing),
                "status": f"error: {str(e)[:80]}",
            })
        finally:
            if tmp_name in adata.obs.columns:
                adata.obs.drop(columns=[tmp_name], inplace=True)

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

    if save_path is not None:
        adata.write(save_path)

    return adata
