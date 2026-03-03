"""
Group-level enrichment testing: ORA and pseudobulk GSEA.

Functions
---------
run_ora                 Over-representation analysis (Fisher exact + BH FDR).
run_gsea_pseudobulk     Pseudobulk GSEA via gseapy prerank.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as stats
from anndata import AnnData
from statsmodels.stats.multitest import multipletests

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _flatten_gene_set_dict(gene_set_dict: dict) -> dict[str, tuple[str, list[str]]]:
    """
    Flatten a possibly two-level gene set dict to {col_name: (category, [genes])}.

    Accepts:
    - Two-level: ``{category: {name: [genes]}}``
    - Flat: ``{name: [genes]}``

    Returns a flat dict keyed by the full path string (e.g. ``"Hallmark/HYPOXIA"``).
    """
    flat: dict[str, tuple[str, list[str]]] = {}
    for k, v in gene_set_dict.items():
        if k == "_meta":
            continue
        if isinstance(v, dict):
            # Two-level
            for name, genes in v.items():
                if isinstance(genes, list):
                    col = f"{k}/{name}"
                    flat[col] = (k, [g for g in genes if isinstance(g, str)])
        elif isinstance(v, list):
            # Flat
            flat[k] = ("", [g for g in v if isinstance(g, str)])
    return flat


def _get_group_genes(adata: AnnData, groupby: str, group: str) -> list[str]:
    """Return genes expressed in at least one cell of the group (count > 0)."""
    mask = adata.obs[groupby] == group
    X = adata.X[mask]
    if sp.issparse(X):
        expressed = np.asarray((X > 0).sum(axis=0)).ravel() > 0
    else:
        expressed = (X > 0).any(axis=0)
    return list(np.array(adata.var_names)[expressed])


# ---------------------------------------------------------------------------
# ORA
# ---------------------------------------------------------------------------


def run_ora(
    adata: AnnData,
    groupby: str,
    gene_set_dict: dict,
    background: list[str] | None = None,
    min_genes: int = 3,
    fdr_method: str = "fdr_bh",
) -> pd.DataFrame:
    """
    Over-representation analysis (ORA) using Fisher exact test.

    For each (group, gene set) pair, tests whether the overlap between the
    gene set and the genes expressed in that group is greater than expected
    by chance, relative to the background gene universe.

    Multiple testing correction is Benjamini-Hochberg (BH) applied across
    all tests within each group.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix. Must have ``adata.obs[groupby]``.
    groupby : str
        Column in ``adata.obs`` defining groups (e.g. ``"leiden"``).
    gene_set_dict : dict
        Gene sets as ``{category: {name: [genes]}}`` or ``{name: [genes]}``.
    background : list[str] or None
        Background gene universe. Defaults to ``adata.var_names``.
    min_genes : int
        Minimum overlap required to include a test (default 3).
    fdr_method : str
        Multiple testing correction method passed to
        ``statsmodels.stats.multitest.multipletests`` (default ``"fdr_bh"``).

    Returns
    -------
    pd.DataFrame
        Long-format results with columns: ``group``, ``category``,
        ``gene_set``, ``n_overlap``, ``n_set``, ``n_background``,
        ``n_group_genes``, ``odds_ratio``, ``p_val``, ``p_adj``.
        Sorted by ``p_adj`` ascending within each group.
    """
    if groupby not in adata.obs.columns:
        raise ValueError(f"groupby={groupby!r} not found in adata.obs.columns")

    bg = list(background) if background is not None else list(adata.var_names)
    bg_set = {g.upper() for g in bg}
    n_background = len(bg_set)

    flat = _flatten_gene_set_dict(gene_set_dict)
    groups = sorted(adata.obs[groupby].unique().tolist())

    rows = []
    for group in groups:
        group_genes = {g.upper() for g in _get_group_genes(adata, groupby, str(group))}
        n_group = len(group_genes)

        p_vals = []
        row_buf = []

        for col_name, (category, genes) in flat.items():
            set_genes = {g.upper() for g in genes if g.upper() in bg_set}
            n_set = len(set_genes)
            if n_set < min_genes:
                continue

            overlap = group_genes & set_genes
            n_overlap = len(overlap)

            # 2x2 contingency:
            # [ overlap,          set_genes - overlap    ]
            # [ group - overlap,  bg - group - set + overlap ]
            a = n_overlap
            b = n_set - n_overlap
            c = n_group - n_overlap
            d = n_background - n_set - n_group + n_overlap
            d = max(d, 0)

            _, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
            odds_ratio = (a * d) / (b * c) if (b * c) > 0 else np.inf

            row_buf.append(
                {
                    "group": group,
                    "category": category,
                    "gene_set": col_name,
                    "n_overlap": n_overlap,
                    "n_set": n_set,
                    "n_background": n_background,
                    "n_group_genes": n_group,
                    "odds_ratio": odds_ratio,
                    "p_val": p,
                }
            )
            p_vals.append(p)

        if not row_buf:
            continue

        # BH correction within group
        _, p_adj, _, _ = multipletests(p_vals, method=fdr_method)
        for row, pa in zip(row_buf, p_adj, strict=False):
            row["p_adj"] = pa
        rows.extend(row_buf)

    if not rows:
        return pd.DataFrame(
            columns=[
                "group",
                "category",
                "gene_set",
                "n_overlap",
                "n_set",
                "n_background",
                "n_group_genes",
                "odds_ratio",
                "p_val",
                "p_adj",
            ]
        )

    result = pd.DataFrame(rows)
    result = result.sort_values(["group", "p_adj"]).reset_index(drop=True)
    return result


# ---------------------------------------------------------------------------
# Pseudobulk GSEA
# ---------------------------------------------------------------------------


def _compute_pseudobulk(adata: AnnData, groupby: str, layer: str | None) -> pd.DataFrame:
    """Aggregate mean expression per group. Returns DataFrame (genes x groups)."""
    groups = sorted(adata.obs[groupby].unique().tolist())
    agg = {}
    for g in groups:
        mask = (adata.obs[groupby] == g).values
        if layer is not None and layer in adata.layers:
            X = adata.layers[layer][mask]
        else:
            X = adata.X[mask]
        if sp.issparse(X):
            X = X.toarray()
        agg[str(g)] = X.mean(axis=0)
    return pd.DataFrame(agg, index=list(adata.var_names))


def _rank_genes(expr: pd.Series, ranking_stat: str, mean_expr: pd.Series) -> pd.Series:
    """Compute a ranking statistic for each gene vs the global mean."""
    if ranking_stat == "logfc":
        # log2 fold change vs mean; add pseudocount
        return np.log2((expr + 1e-9) / (mean_expr + 1e-9))
    elif ranking_stat == "zscore":
        std = mean_expr.std()
        if std == 0:
            return pd.Series(0.0, index=expr.index)
        return (expr - mean_expr) / std
    elif ranking_stat == "tstat":
        # Approximate t-stat: (group - mean) / std(mean)
        std = mean_expr.std()
        if std == 0:
            return pd.Series(0.0, index=expr.index)
        return (expr - mean_expr) / std
    else:
        raise ValueError(f"ranking_stat={ranking_stat!r}. Choose from 'logfc', 'tstat', 'zscore'.")


def run_gsea_pseudobulk(
    adata: AnnData,
    groupby: str,
    gene_set_dict: dict,
    layer: str | None = None,
    method: str = "prerank",
    ranking_stat: str = "logfc",
    min_genes: int = 10,
    n_permutations: int = 1000,
) -> pd.DataFrame:
    """
    Pseudobulk GSEA: aggregate expression per group, rank genes, run fgsea/prerank.

    Requires ``gseapy`` (optional dependency). Install with::

        pip install gseapy

    Steps:

    1. Aggregate expression: mean per group → pseudobulk matrix (genes x groups).
    2. Compute ranking statistic (log2FC, z-score, or t-stat) relative to the
       mean expression across all groups.
    3. Pass ranked gene list to ``gseapy.prerank`` for each group.
    4. Collect results; apply Benjamini-Hochberg FDR across all gene sets per group.

    Parameters
    ----------
    adata : AnnData
        Annotated data. Must have ``adata.obs[groupby]``.
    groupby : str
        Column in ``adata.obs`` defining groups.
    gene_set_dict : dict
        Gene sets as ``{category: {name: [genes]}}`` or ``{name: [genes]}``.
    layer : str or None
        Layer to use for expression. Defaults to ``adata.X``.
    method : str
        GSEA method: ``"prerank"`` (default) or ``"fgsea"`` (also via gseapy).
    ranking_stat : str
        Statistic used to rank genes: ``"logfc"`` (default), ``"zscore"``, ``"tstat"``.
    min_genes : int
        Minimum genes in a set after intersecting with ranked genes (default 10).
    n_permutations : int
        Number of permutations for enrichment scoring (default 1000).

    Returns
    -------
    pd.DataFrame
        Long-format results with columns: ``group``, ``category``, ``gene_set``,
        ``NES``, ``p_val``, ``p_adj``, ``lead_edge_genes``.
        Sorted by ``p_adj`` ascending within each group.
    """
    try:
        import gseapy
    except ImportError as exc:
        raise ImportError(
            "gseapy is required for run_gsea_pseudobulk. Install it with: pip install gseapy"
        ) from exc

    if groupby not in adata.obs.columns:
        raise ValueError(f"groupby={groupby!r} not found in adata.obs.columns")

    # Build flat gene set dict for gseapy (name -> [genes])
    flat = _flatten_gene_set_dict(gene_set_dict)
    gseapy_gene_sets = {col: genes for col, (_, genes) in flat.items()}

    # Pseudobulk expression
    pseudobulk = _compute_pseudobulk(adata, groupby, layer)
    mean_expr = pseudobulk.mean(axis=1)

    groups = sorted(pseudobulk.columns.tolist())
    rows = []

    for group in groups:
        expr = pseudobulk[group]
        ranked = _rank_genes(expr, ranking_stat, mean_expr).dropna().sort_values(ascending=False)

        if ranked.empty:
            continue

        try:
            if method in ("prerank", "fgsea"):
                gs_res = gseapy.prerank(
                    rnk=ranked,
                    gene_sets=gseapy_gene_sets,
                    min_size=min_genes,
                    permutation_num=n_permutations,
                    outdir=None,
                    no_plot=True,
                    verbose=False,
                )
                res_df = gs_res.res2d if hasattr(gs_res, "res2d") else pd.DataFrame()
            else:
                raise ValueError(f"method={method!r}. Choose 'prerank' or 'fgsea'.")
        except Exception as e:
            rows.append(
                {
                    "group": group,
                    "category": "",
                    "gene_set": "ERROR",
                    "NES": np.nan,
                    "p_val": np.nan,
                    "p_adj": np.nan,
                    "lead_edge_genes": str(e)[:120],
                }
            )
            continue

        if res_df.empty:
            continue

        p_vals = []
        row_buf = []
        for col_name, (category, _) in flat.items():
            if col_name not in res_df.index:
                continue
            row_entry = res_df.loc[col_name]
            nes = float(row_entry.get("NES", np.nan))
            p_val = float(row_entry.get("NOM p-val", row_entry.get("pval", np.nan)))
            lead = row_entry.get("Lead_genes", row_entry.get("lead_genes", ""))
            row_buf.append(
                {
                    "group": group,
                    "category": category,
                    "gene_set": col_name,
                    "NES": nes,
                    "p_val": p_val,
                    "lead_edge_genes": lead
                    if isinstance(lead, str)
                    else ";".join(lead)
                    if lead
                    else "",
                }
            )
            p_vals.append(p_val if not np.isnan(p_val) else 1.0)

        if not row_buf:
            continue

        _, p_adj, _, _ = multipletests(p_vals, method="fdr_bh")
        for row_entry, pa in zip(row_buf, p_adj, strict=False):
            row_entry["p_adj"] = pa
        rows.extend(row_buf)

    if not rows:
        return pd.DataFrame(
            columns=["group", "category", "gene_set", "NES", "p_val", "p_adj", "lead_edge_genes"]
        )

    result = pd.DataFrame(rows)
    result = result.sort_values(["group", "p_adj"]).reset_index(drop=True)
    return result
