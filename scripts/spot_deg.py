# -*- coding: utf-8 -*-
"""
DEG per sample (1 vs rest within user-defined sets) using a normalized+log1p view of adata.raw.X.
- Keeps ALL genes (no HVG filter)
- Top 5 significant genes per sample (padj <= PADJ_MAX and |logFC| >= MIN_ABS_LOGFC)
- Dotplot per compartment (epithelial, stromal, immune)
"""

from pathlib import Path
from typing import Dict, List, Sequence, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import warnings
from scipy import sparse

# -----------------------------
# Config
# -----------------------------
ADATA_PATH = Path("results/scvi.leiden.phenotyped.h5ad")
cluster_key: str = "cluster"
library_id_key: str = "library_id"
library_ids: Dict[str, List[str]] = {
    "P1_NAT": ["PT01-1_NAT", "PT01-3_NAT", "PT01-5_NAT"],
    "P1_TUM": ["PT01-2_TUM", "PT01-4_TUM", "PT01-6_TUM"],
    "P3_NAT": ["P3_S1_BL_NT_1_iF10_S7"],
    "P3_TUM": ["P3_S2_BL_T_2_iG10_S8"],
    "P5_NAT": ["Pat5_Samp1_BL_2_iD10_S5", "Pat_5_BL_Tu_Samp2_A1_iD9_S4"],
    "P1_LN": ["Pat_1_LN_PRT_Samp8_D1_iC9_S3"],
    "P3_LN": ["Pat_3_LN_PRT_Samp_8_D1_iB9_S2", "Pat_3_LN_Samp_7_D1_iA9_S1"],
}

# DEG parameters
DE_METHOD = "wilcoxon"          # robust and common with log1p data
PADJ_MAX = 0.05
MIN_ABS_LOGFC = 0.5             # on log1p scale; adjust if needed
N_TOP = 5

# Plot parameters
DOT_CMAP = "RdBu_r"
FIGSIZE = (12, 8)
DOTSIZE_MIN, DOTSIZE_MAX = 40, 300
OUT_DIR = Path("figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# -----------------------------
# Utilities
# -----------------------------
def _compartment(label: str) -> str:
    l = label.lower()
    if any(k in l for k in [
        "enterocyte", "goblet", "enteroendocrine", "crypt stem", "regenerative",
        "neuroendocrine", "metabolic / stress", "weak epithelial", "epithelium", "tumor"
    ]):
        return "epithelial"
    if any(k in l for k in ["fibroblast", "smooth muscle", "pericyte", "ecm mix", "matrix remodeling", "ecm high"]):
        return "stromal"
    if any(k in l for k in ["b cell", "plasma", "gc", "myeloid", "mast", "nk", "t cell", "dendritic"]):
        return "immune"
    return "unknown"

def _validate(adata: sc.AnnData) -> None:
    assert cluster_key in adata.obs, f"obs lacks {cluster_key}"
    assert library_id_key in adata.obs, f"obs lacks {library_id_key}"
    assert adata.raw is not None and adata.raw.shape[1] > 0, "adata.raw is required and must contain genes"
    if not pd.api.types.is_categorical_dtype(adata.obs[cluster_key]):
        adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")
    adata.obs[library_id_key] = adata.obs[library_id_key].astype(str)

def _mask_compartment(adata: sc.AnnData, comp: str) -> np.ndarray:
    return adata.obs[cluster_key].astype(str).map(_compartment).values == comp

def _subset_present_samples(adata: sc.AnnData, desired: Sequence[str]) -> List[str]:
    present = set(adata.obs[library_id_key].unique().tolist())
    return [s for s in desired if s in present]

def _normalized_log1p_view_from_raw(A: sc.AnnData) -> sc.AnnData:
    """
    Build a temporary AnnData whose X is adata.raw.X normalized and log1p-transformed.
    Does not mutate the input object.
    """
    assert A.raw is not None
    Xraw = A.raw.X
    # Create a lightweight AnnData with raw counts in X
    B = sc.AnnData(
        X=Xraw.copy() if isinstance(Xraw, np.ndarray) else Xraw.copy(),
        obs=A.obs.copy(),
        var=A.raw.var.copy(),
    )
    # library-size normalize and log1p
    sc.pp.normalize_total(B, target_sum=1e4, inplace=True)
    sc.pp.log1p(B)  # natural log(1 + normalized)
    return B

def _rank_genes_1vRest_normlog(
    As: sc.AnnData,
    groupby: str,
    target_group: str,
) -> pd.DataFrame:
    """
    Run 1 vs rest DE on a normalized+log1p view of As.raw.X.
    Returns a dataframe for the target group with logFC on log1p scale.
    """
    B = _normalized_log1p_view_from_raw(As)
    if not pd.api.types.is_categorical_dtype(B.obs[groupby]):
        B.obs[groupby] = B.obs[groupby].astype("category")
    sc.tl.rank_genes_groups(
        B,
        groupby=groupby,
        groups=[target_group],
        reference="rest",
        method=DE_METHOD,
        use_raw=False,             # we already put norm+log in X
        n_genes=B.shape[1],        # all genes
    )
    df = sc.get.rank_genes_groups_df(B, group=target_group)
    keep = ["names", "logfoldchanges", "scores", "pvals_adj"]
    df = df[[c for c in keep if c in df.columns]]
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["logfoldchanges", "pvals_adj"])
    return df

def _percent_expressing_raw(A: sc.AnnData, sample_name: str, genes: Sequence[str]) -> pd.Series:
    """
    Compute percent expressing from A.raw.X for one sample (within the current compartment).
    """
    assert A.raw is not None
    mask_cells = (A.obs[library_id_key].astype(str) == sample_name).values
    if mask_cells.sum() == 0:
        return pd.Series(index=genes, dtype=float)
    Xs = A.raw[mask_cells, :].X
    if sparse.issparse(Xs):
        nz = (Xs > 0).astype(float)
        pct = 100.0 * (nz.mean(axis=0).A1)
    else:
        pct = 100.0 * (Xs > 0).mean(axis=0)
    s = pd.Series(pct, index=A.raw.var_names)
    return s.reindex(genes)

def _collect_deg_for_compartment(
    adata: sc.AnnData,
    comp: str,
    library_sets: Dict[str, List[str]],
    n_top: int = N_TOP,
    padj_max: float = PADJ_MAX,
    min_abs_logfc: float = MIN_ABS_LOGFC,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    For a compartment:
      - For each user set with >= 2 samples, run 1vRest per sample on a norm+log view of raw
      - Keep up to n_top genes per sample that satisfy padj_max and min_abs_logfc
      - Return df_lfc (logFC) and df_pct (% expressing from raw) with NaNs elsewhere
    """
    mask = _mask_compartment(adata, comp)
    if mask.sum() == 0:
        return pd.DataFrame(), pd.DataFrame()

    A = adata[mask].copy()
    A.obs[library_id_key] = A.obs[library_id_key].astype(str).astype("category")

    sample2genes: Dict[str, List[str]] = {}
    lfc_store: Dict[Tuple[str, str], float] = {}

    for set_name, samples in library_sets.items():
        samples_present = _subset_present_samples(A, samples)
        if len(samples_present) < 2:
            continue
        As = A[A.obs[library_id_key].isin(samples_present)].copy()
        As.obs[library_id_key] = As.obs[library_id_key].astype("category")

        for samp in samples_present:
            df = _rank_genes_1vRest_normlog(As, groupby=library_id_key, target_group=samp)
            sig = df[(df["pvals_adj"] <= padj_max) & (df["logfoldchanges"].abs() >= min_abs_logfc)].copy()
            if sig.empty:
                continue
            # Remove unwanted genes before selecting top N
            bad_prefixes = ("RPS", "RPL", "MT-", "MT_", "HB", "HBA", "HBB", "HBD", "HBG", "HBM", "HBQ", "HBZ")
            sig = sig[~sig["names"].str.upper().str.startswith(bad_prefixes)]

            # Now pick top N from filtered genes
            sig["abs_lfc"] = sig["logfoldchanges"].abs()
            sig = sig.sort_values("abs_lfc", ascending=False).head(n_top)
            genes = sig["names"].astype(str).tolist()

            sample2genes[samp] = genes
            for g, lfc in zip(sig["names"], sig["logfoldchanges"]):
                lfc_store[(str(g), samp)] = float(lfc)
    
    if not sample2genes:
        return pd.DataFrame(), pd.DataFrame()

    genes = sorted({g for gl in sample2genes.values() for g in gl})
    samples_all = sorted(sample2genes.keys())

    df_lfc = pd.DataFrame(np.nan, index=genes, columns=samples_all, dtype=float)
    df_pct = pd.DataFrame(np.nan, index=genes, columns=samples_all, dtype=float)

    for s in samples_all:
        genes_s = sample2genes[s]
        pct_s = _percent_expressing_raw(A, s, genes_s)
        for g in genes_s:
            if (g, s) in lfc_store:
                df_lfc.at[g, s] = lfc_store[(g, s)]
            if g in pct_s.index and pd.notnull(pct_s.loc[g]):
                df_pct.at[g, s] = float(pct_s.loc[g])

    # make dots readable
    df_pct = df_pct.clip(lower=5.0, upper=100.0)

    # --- NEW BLOCK 1: remove ribosomal, mitochondrial, and hemoglobin genes ---
    bad_prefixes = ("RPS", "RPL", "MT-", "MT_", "HB", "HBA", "HBB", "HBD", "HBG", "HBM", "HBQ", "HBZ")
    mask_valid = ~df_lfc.index.str.upper().str.startswith(bad_prefixes)
    df_lfc = df_lfc.loc[mask_valid]
    df_pct = df_pct.loc[mask_valid]

    # --- NEW BLOCK 2: reorder genes by their strongest (abs) logFC sample ---
    dominant_sample = df_lfc.abs().idxmax(axis=1)
    order = []
    for s in samples_all:
        genes_s = dominant_sample[dominant_sample == s].index.tolist()
        genes_s = sorted(genes_s, key=lambda g: abs(df_lfc.loc[g, s]), reverse=True)
        order.extend(genes_s)
    order = [g for g in order if g in df_lfc.index]

    df_lfc = df_lfc.loc[order]
    df_pct = df_pct.loc[order]

    return df_lfc, df_pct


def _dotplot_from_tables(
    df_lfc: pd.DataFrame,
    df_pct: pd.DataFrame,
    title: str,
    out_png: Path,
    cmap: str = DOT_CMAP,
    figsize: Tuple[int, int] = FIGSIZE,
    dot_min: float = DOTSIZE_MIN,
    dot_max: float = DOTSIZE_MAX,
) -> None:
    """Dotplot of logFC and % expressing; draw only non-NaN entries."""
    if df_lfc.empty:
        warnings.warn(f"No DE results to plot for: {title}")
        return
    df_pct = df_pct.loc[df_lfc.index, df_lfc.columns]

    V, S = df_lfc.values, df_pct.values
    mask = ~np.isnan(V) & ~np.isnan(S)
    if not mask.any():
        warnings.warn(f"No valid entries to plot for: {title}")
        return

    sizes = dot_min + (np.clip(S, 0, 100) / 100.0) * (dot_max - dot_min)
    absvals = np.abs(V[mask])
    vmax = np.nanpercentile(absvals, 98)
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = np.nanmax(absvals) if np.isfinite(np.nanmax(absvals)) else 1.0
    vmin = -vmax

    xs, ys = np.meshgrid(np.arange(df_lfc.shape[1]), np.arange(df_lfc.shape[0]))
    xs, ys = xs[mask], ys[mask]

    fig, ax = plt.subplots(figsize=figsize)
    sca = ax.scatter(xs, ys, c=V[mask], s=sizes[mask], cmap=cmap, vmin=vmin, vmax=vmax,
                     edgecolors="none", alpha=0.95)
    ax.set_xticks(np.arange(df_lfc.shape[1]))
    ax.set_xticklabels(df_lfc.columns.tolist(), rotation=90)
    ax.set_yticks(np.arange(df_lfc.shape[0]))
    ax.set_yticklabels(df_lfc.index.tolist())
    ax.set_title(title)
    ax.invert_yaxis()
    ax.set_xlim(-0.5, df_lfc.shape[1] - 0.5)
    ax.set_ylim(-0.5, df_lfc.shape[0] - 0.5)

    cbar = fig.colorbar(sca, ax=ax, shrink=0.85)
    cbar.set_label("log fold change (1 vs rest; norm+log1p)")

    for p in [25, 50, 75, 100]:
        plt.scatter([], [], s=dot_min + (p/100)*(dot_max - dot_min), c="k", label=f"{p}%")
    ax.legend(scatterpoints=1, frameon=False, title="% expressing",
              loc="upper right", bbox_to_anchor=(1.22, 1.0))

    for i, sample in enumerate(df_lfc.columns):
        if i % 2 == 0:
            ax.axvspan(i - 0.5, i + 0.5, color='lightgrey', alpha=0.1, zorder=0)

    fig.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

# -----------------------------
# Main
# -----------------------------
def main():
    adata = sc.read(ADATA_PATH)
    _validate(adata)

    jobs = [
        ("epithelial", OUT_DIR / "deg_dotplot_epithelial.png"),
        ("stromal",    OUT_DIR / "deg_dotplot_stromal.png"),
        ("immune",     OUT_DIR / "deg_dotplot_immune.png"),
    ]

    for comp, out_png in jobs:
        print(f"[{comp}] DE on normalized+log1p view of adata.raw.X; selecting up to {N_TOP} significant genes per sample")
        df_lfc, df_pct = _collect_deg_for_compartment(
            adata=adata,
            comp=comp,
            library_sets=library_ids,
            n_top=N_TOP,
            padj_max=PADJ_MAX,
            min_abs_logfc=MIN_ABS_LOGFC,
        )
        if df_lfc.empty:
            print(f"[{comp}] no significant genes found (given thresholds) or no sets with >=2 samples. Skipping.")
            continue

        base = out_png.with_suffix("")
        df_lfc.to_csv(base.as_posix() + "_lfc.csv")
        df_pct.to_csv(base.as_posix() + "_pct.csv")

        _dotplot_from_tables(
            df_lfc=df_lfc,
            df_pct=df_pct,
            title=f"{comp.capitalize()} – top {N_TOP} DE genes per sample (1 vs rest; norm+log1p on raw)",
            out_png=out_png,
            cmap=DOT_CMAP,
            figsize=FIGSIZE,
        )
        print(f"[{comp}] saved {out_png}")

if __name__ == "__main__":
    main()
