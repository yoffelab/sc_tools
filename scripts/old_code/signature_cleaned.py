# ============================
# Minimal, modular signatures
# ============================
import scanpy as sc
import json

sc.settings.set_figure_params(dpi=300, dpi_save = 400) 

# Assuming 'data.json' contains valid JSON data
with open('metadata/gene_signatures.json', 'r') as file:
    spatial_signatures = json.load(file)
print(spatial_signatures)


from __future__ import annotations
from typing import Mapping, Any, Iterable, Tuple, List, Dict
import numpy as np, pandas as pd, scipy.sparse as sp
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu


# ---------- utilities ----------
def _flatten_nested_dict(d: Mapping[str, Any], prefix: Tuple[str,...]=()) -> List[Tuple[Tuple[str,...], List[str]]]:
    out = []
    for k, v in d.items():
        if isinstance(v, Mapping):
            out.extend(_flatten_nested_dict(v, prefix + (k,)))
        elif isinstance(v, list):
            out.append((prefix + (k,), [g for g in v if isinstance(g, str)]))
    return out

def _index_genes(genes: Iterable[str], var_names: np.ndarray) -> Tuple[np.ndarray, List[str], List[str]]:
    vm = {g.upper(): i for i, g in enumerate(var_names)}
    idx, present, missing = [], [], []
    for g in genes:
        key = g.upper()
        if key in vm:
            idx.append(vm[key]); present.append(var_names[vm[key]])
        else:
            missing.append(g)
    return np.array(sorted(set(idx))), present, missing

def _get_matrix(adata, layer: str|None, use_raw: bool, log1p: bool, clip_min: float|None):
    assert not (use_raw and (layer is not None)), "Choose use_raw or layer, not both"
    if use_raw:
        assert adata.raw is not None, "adata.raw missing but use_raw=True"
        X, varn = adata.raw.X, np.array(adata.raw.var_names, str)
    elif layer is not None:
        assert layer in adata.layers, f"Layer {layer} not in adata.layers"
        X, varn = adata.layers[layer], np.array(adata.var_names, str)
    else:
        X, varn = adata.X, np.array(adata.var_names, str)
    if sp.issparse(X): X = X.tocsr()
    else: X = np.asarray(X, float)
    if log1p:
        if sp.issparse(X):
            if clip_min is not None:
                m = X.data < clip_min
                if np.any(m): X.data[m] = clip_min
            else:
                assert np.all(X.data >= 0), "Negative values with log1p"
            X.data = np.log1p(X.data)
        else:
            if clip_min is not None: X = np.clip(X, clip_min, None)
            else: assert np.all(X >= 0), "Negative values with log1p"
            X = np.log1p(X)
    return X, varn

def _zscore_cols(X):
    eps = 1e-8
    if sp.issparse(X):
        n = X.shape[0]
        s = np.asarray(X.sum(axis=0)).ravel()
        ss = np.asarray(X.power(2).sum(axis=0)).ravel()
        mu = s / max(n, 1)
        var = np.maximum(ss / max(n, 1) - mu**2, 0.0)
        std = np.sqrt(var) + eps
        Xc = X.tocsc(copy=True)
        for j in range(Xc.shape[1]):
            a, b = Xc.indptr[j], Xc.indptr[j+1]
            if b > a: Xc.data[a:b] = (Xc.data[a:b] - mu[j]) / std[j]
        A = Xc.toarray(); A += (-mu / std)
        return A
    mu = X.mean(axis=0, dtype=float)
    sd = X.std(axis=0, ddof=0, dtype=float) + eps
    return (X - mu) / sd

def score_signatures_nested(
    adata,
    signatures_nested: Dict[str, Any],
    prefix: str = "sig",
    use_raw: bool = True,
    layer: str | None = None,
    log1p: bool = True,
    clip_min: float | None = 0.0,
    zscore_per_gene: bool = True,
    min_genes: int = 3,
) -> pd.DataFrame:
    X, var_names = _get_matrix(adata, layer, use_raw, log1p, clip_min)
    leaves = _flatten_nested_dict(signatures_nested)
    assert len(leaves) > 0, "No gene lists in signatures_nested"
    results, report = {}, []
    for path, genes in leaves:
        idx, present, missing = _index_genes(genes, var_names)
        col = f"{prefix}:{'/'.join(path)}"
        if idx.size < min_genes:
            results[col] = np.full(adata.n_obs, np.nan)
            report.append({"column": col, "n_present": int(idx.size), "status": "skipped"})
            continue
        X_sub = X[:, idx]
        if zscore_per_gene: X_sub = _zscore_cols(X_sub)
        score = np.asarray(X_sub.mean(axis=1)).ravel()
        results[col] = score
        report.append({"column": col, "n_present": int(idx.size), "status": "ok"})
    df = pd.DataFrame(results, index=adata.obs_names)
    for c in df.columns: adata.obs[c] = df[c].astype(float)
    adata.uns[prefix] = {"report": pd.DataFrame(report)}
    return df

# ---------- simple plotting ----------
def present_cols(adata, cols):
    return [c for c in cols if c in adata.obs.columns]

def plot_panel(
    adata_sub,
    cols,
    library_id,                # <-- add this
    title=None,
    cmap="coolwarm",
    vmin=-2, vmax=2, vcenter=0,
    ncols=4,
):
    cols = [c for c in cols if c in adata_sub.obs.columns]
    if len(cols) == 0:
        print(f"Warning: No columns to plot for {library_id}, skipping...")
        return  # Skip instead of failing
    assert len(cols) > 0, "No columns to plot"
    # critical: specify the library
    sc.pl.spatial(
        adata_sub,
        color=cols,
        library_id=library_id,   # <-- tell Scanpy which image to use
        use_raw=False,
        ncols=ncols,
        cmap=cmap,
        vmin=vmin, vmax=vmax, vcenter=vcenter,
        frameon=False,
        show=False,
    )
    if title:
        import matplotlib.pyplot as plt
        plt.suptitle(title)

def compare_core_tls(df, solidity_order=["Normal", "Non-Solid", "Solid"]):
    """
    Performs two-sided Mann-Whitney U tests comparing Core vs TLS
    for each signature column (columns starting with 'sig:/'),
    within each solidity group.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain columns:
            - 'solidity_type'
            - 'architecture_type'
            - signature columns starting with 'sig:/'
    solidity_order : list or None
        Order of solidity groups (e.g., ["Normal", "Non-Solid", "Solid"])

    Returns
    -------
    pd.DataFrame
        index: (solidity_type, signature)
        columns: ['pval', 'n_core', 'n_tls']
    """

    # Identify signature columns
    sig_cols = [c for c in df.columns if c.startswith("sig:")]
    assert len(sig_cols) > 0, "No signature columns starting with 'sig:' found."

    if solidity_order is None:
        solidity_order = df["solidity_type"].unique().tolist()

    results = []

    for sol in solidity_order:
        df_sol = df[df["solidity_type"] == sol]

        df_core = df_sol[df_sol["architecture_type"] == "Core"]
        df_tls  = df_sol[df_sol["architecture_type"] == "TLS"]

        # Skip groups with insufficient data
        if len(df_core) < 3 or len(df_tls) < 3:
            for sig in sig_cols:
                results.append([sol, sig, np.nan, len(df_core), len(df_tls)])
            continue

        for sig in sig_cols:
            x = df_core[sig].dropna()
            y = df_tls[sig].dropna()

            if len(x) < 3 or len(y) < 3:
                pval = np.nan
            else:
                _, pval = mannwhitneyu(x, y, alternative="two-sided")

            results.append([sol, sig, pval, len(x), len(y)])

    res_df = pd.DataFrame(results, columns=["solidity_type", "signature", "pval", "n_core", "n_tls"])
    res_df.set_index(["solidity_type", "signature"], inplace=True)

    return res_df
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def plot_signature_boxplots(df, stats_df, solidity_order=None, figsize=(6, 4)):

    # Extract signature columns (note: using "sig:" not "sig:/")
    sig_cols = [c for c in df.columns if c.startswith("sig:")]
    assert len(sig_cols) > 0, "No signature columns found (columns starting with 'sig:')."

    if solidity_order is None:
        solidity_order = df["solidity_type"].unique().tolist()

    for sol in solidity_order:

        df_sol_full = df[df["solidity_type"] == sol]

        fig, axes = plt.subplots(
            nrows=1,
            ncols=len(sig_cols),
            figsize=(figsize[0] * len(sig_cols), figsize[1]),
        )

        if len(sig_cols) == 1:
            axes = [axes]

        for ax, sig in zip(axes, sig_cols):

            # Filter to Core/TLS only
            df_sol = df_sol_full[df_sol_full["architecture_type"].isin(["Core", "TLS"])].copy()

            cats = sorted(df_sol["architecture_type"].unique().tolist())

            # -------------------------
            # CASE 1: No Core or TLS present
            # -------------------------
            if len(cats) == 0:
                ax.text(0.5, 0.5, f"{sol}\n{sig}\nNo Core or TLS", ha="center")
                ax.set_axis_off()
                continue

            # -------------------------
            # CASE 2: Only one category present
            # -------------------------
            if len(cats) == 1:
                single = cats[0]
                sns.boxplot(
                    data=df_sol,
                    x="architecture_type",
                    y=sig,
                    order=[single],       # ONLY one category plotted
                    color="#a6cee3",
                    ax=ax
                )
                sns.stripplot(
                    data=df_sol,
                    x="architecture_type",
                    y=sig,
                    order=[single],
                    color="black",
                    size=2,
                    ax=ax
                )

                ax.text(
                    0.5, 1.05,
                    "Only one group present",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                    fontsize=9,
                )
                ax.set_title(f"{sol} — {sig}", fontsize=10)
                continue

            # -------------------------
            # CASE 3: Both Core AND TLS exist → full statistical test
            # -------------------------
            sns.boxplot(
                data=df_sol,
                x="architecture_type",
                y=sig,
                order=["Core", "TLS"],
                color="#a6cee3",
                ax=ax
            )

            sns.stripplot(
                data=df_sol,
                x="architecture_type",
                y=sig,
                order=["Core", "TLS"],
                color="black",
                size=2,
                ax=ax
            )

            # Add p-value
            try:
                pval = stats_df.loc[(sol, sig), "pval"]
            except KeyError:
                pval = np.nan

            if np.isnan(pval):
                label = "n.s."
            elif pval < 0.05:
                label = f"p = {pval:.3g}"
            else:
                label = "n.s."

            ax.text(
                0.5, 1.05,
                label,
                ha="center",
                transform=ax.transAxes,
                fontsize=10,
                fontweight="bold"
            )

            ax.set_title(f"{sol} — {sig}", fontsize=10)

        plt.tight_layout()
        plt.show()

adata = sc.read('results/adata.annotation.masked.h5ad')
# 1) Score once
scores = score_signatures_nested(
    adata,
    signatures_nested=spatial_signatures,
    prefix="sig",
    use_raw=True,   # or layer="lognorm"
    log1p=True,
    zscore_per_gene=True,
)

# 2) Z-score each signature across spots for consistent contrast
for c in [c for c in adata.obs.columns if c.startswith("sig:")]:
    x = adata.obs[c].values
    adata.obs[c + "_z"] = (x - np.nanmean(x)) / (np.nanstd(x) + 1e-8)

# 3) Plot per slide with concise panels
myeloid_panel = [
    "sig:Myeloid/Macrophage_M1_z",
    "sig:Myeloid/Macrophage_M2_z",
    "sig:Myeloid/Macrophage_Repair_z",
    "sig:Myeloid/Dendritic_cDC1_z",
    "sig:Myeloid/Dendritic_cDC2_z",
    "sig:Myeloid/Neutrophil_Core_z",
]
tumor_panel = [
    "sig:CRC_Tumor/Epithelial_Core_z",
    "sig:CRC_Tumor/Stem_WNT_MYC_z",
    "sig:CRC_Tumor/TGFb_EMT_Invasion_z",
    "sig:CRC_Tumor/MSI_Immune_HLAII_z",
    "sig:CRC_Tumor/CIN_Mitotic_z",
    "sig:Proliferation/Cell_Cycle_z",
]

lymphoid_panel = [
    "sig:Lymphoid/CD4_Naive_z",
    "sig:Lymphoid/CD4_Th1_z",
    "sig:Lymphoid/CD4_Th2_z",
    "sig:Lymphoid/CD4_Th17_Gut_z",
    "sig:Lymphoid/CD4_Treg_z",
    "sig:Lymphoid/CD4_Exhausted_z",
    "sig:Lymphoid/CD8_Naive_z",
    "sig:Lymphoid/CD8_Effector_z",
    "sig:Lymphoid/CD8_Exhausted_z",
    "sig:Lymphoid/CD8_TRM_z",
    "sig:Lymphoid/CD8_Effector_KIR_z",
    "sig:Lymphoid/NK_Cytotoxic_z",
    "sig:Lymphoid/NK_Regulatory_z",
    "sig:Lymphoid/NKT_z",
    "sig:Lymphoid/ILC2_Repair_z",
    "sig:Lymphoid/ILC3_Mucosal_z",
    "sig:Lymphoid/B_Core_z",
    "sig:Lymphoid/B_Memory_z",
    "sig:Lymphoid/B_Naive_z",
    "sig:Lymphoid/B_GerminalCenter_z",
    "sig:Lymphoid/Plasma_Core_z",
    "sig:Lymphoid/Plasma_IgA_z",
    "sig:Lymphoid/Plasma_IgG_z",
    "sig:Lymphoid/Plasma_IgM_IgD_z",
    "sig:Lymphoid/Plasmablast_z",
]

process_panel = [
    "sig:CRC_Tumor/Epithelial_Core_z",
    "sig:Proliferation/Cell_Cycle_z",
    "sig:Processes/Angiogenesis_z",
    "sig:Processes/ECM_Remodeling_z",
    "sig:Myeloid/Macrophage_M1_z",
    "sig:Myeloid/Macrophage_M2_z",
    "sig:Myeloid/Macrophage_Repair_z",
]

from tqdm import tqdm
for lib in tqdm(adata.obs["library_id"].unique()):
    a = adata[adata.obs["library_id"] == lib].copy()
    plot_panel(a, process_panel, library_id=lib, title=f"{lib} Processes", ncols = 7)
    plt.subplots_adjust(
        left=0.05,   # reduce whitespace on left
        right=0.95,  # reduce whitespace on right
        top=0.90,    # leave room for suptitle
        bottom=0.05, # reduce whitespace at bottom
        wspace=0.1,  # horizontal spacing between subplots
        hspace=0.1   # vertical spacing between subplots
    )
    plt.savefig(f"figures/spatial/{lib}_processes.png", bbox_inches="tight"); plt.close()

for lib in tqdm(adata.obs["library_id"].unique()):
    a = adata[adata.obs["library_id"] == lib].copy()
    plot_panel(a, lymphoid_panel, library_id=lib, title=f"{lib} Processes", ncols = 7)
    plt.subplots_adjust(
        left=0.05,   # reduce whitespace on left
        right=0.95,  # reduce whitespace on right
        top=0.90,    # leave room for suptitle
        bottom=0.05, # reduce whitespace at bottom
        wspace=0.1,  # horizontal spacing between subplots
        hspace=0.1   # vertical spacing between subplots
    )
    plt.savefig(f"figures/spatial/{lib}_lymphoid.png", bbox_inches="tight"); plt.close()




# Statistical Comparison
# 1. Run statistics
df = adata.obs.copy()

# restore or calculate signature columns
# now filter after signatures exist
df = df[df['architecture_type'].isin(['Core','TLS'])].copy()
df['architecture_type'] = pd.Categorical(df['architecture_type'], categories=['Core','TLS'])

# Now compute stats
stats_df = compare_core_tls(df, solidity_order=["Normal", "Non-Solid", "Solid"])

# Now plot
plot_signature_boxplots(df, stats_df, solidity_order=["Normal", "Non-Solid", "Solid"])






import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu


def compare_core_tls_by_solidity(df, sig_prefix="sig:", sig_suffix="_z", solidity_order=None, min_samples=3):
    """
    Performs two-sided Mann-Whitney U tests comparing Core vs TLS
    for each signature column, within each solidity group.
    ALSO performs pairwise comparisons between solidity types within each architecture.
    
    Parameters
    ----------
    df : pd.DataFrame
        Must contain columns:
            - 'solidity_type'
            - 'architecture_type'
            - signature columns starting with sig_prefix
    sig_prefix : str
        Prefix for signature columns (default: "sig:")
    solidity_order : list or None
        Order of solidity groups (e.g., ["Normal", "Non-Solid", "Solid"])
    min_samples : int
        Minimum number of samples required in each group for comparison
    
    Returns
    -------
    tuple: (core_vs_tls_df, pairwise_solidity_df)
        core_vs_tls_df: Core vs TLS comparisons within each solidity
        pairwise_solidity_df: Pairwise solidity comparisons within each architecture
    """
    # Identify signature columns
    sig_cols = [c for c in df.columns if c.startswith(sig_prefix) and c.endswith(sig_suffix)]
    
    if len(sig_cols) == 0:
        raise ValueError(f"No signature columns found starting with '{sig_prefix}'")
    
    print(f"Found {len(sig_cols)} signature columns")
    
    # Determine solidity order
    if solidity_order is None:
        solidity_order = sorted(df["solidity_type"].dropna().unique())
    
    # ===== Part 1: Core vs TLS within each solidity =====
    results_core_tls = []
    
    for sol in solidity_order:
        df_sol = df[df["solidity_type"] == sol].copy()
        df_core = df_sol[df_sol["architecture_type"] == "Core"]
        df_tls = df_sol[df_sol["architecture_type"] == "TLS"]
        
        n_core_total = len(df_core)
        n_tls_total = len(df_tls)
        
        if n_core_total == 0 and n_tls_total == 0:
            status = "no_data"
        elif n_core_total == 0:
            status = "no_core"
        elif n_tls_total == 0:
            status = "no_tls"
        else:
            status = "ok"
        
        for sig in sig_cols:
            x_core = df_core[sig].dropna()
            x_tls = df_tls[sig].dropna()
            
            n_core = len(x_core)
            n_tls = len(x_tls)
            
            mean_core = x_core.mean() if n_core > 0 else np.nan
            mean_tls = x_tls.mean() if n_tls > 0 else np.nan
            
            if n_core >= min_samples and n_tls >= min_samples:
                try:
                    _, pval = mannwhitneyu(x_core, x_tls, alternative="two-sided")
                    test_status = status
                except Exception as e:
                    pval = np.nan
                    test_status = f"error: {str(e)}"
            else:
                pval = np.nan
                if n_core < min_samples and n_tls < min_samples:
                    test_status = "insufficient_both"
                elif n_core < min_samples:
                    test_status = "insufficient_core"
                else:
                    test_status = "insufficient_tls"
            
            results_core_tls.append({
                "solidity_type": sol,
                "signature": sig,
                "pval": pval,
                "n_core": n_core,
                "n_tls": n_tls,
                "mean_core": mean_core,
                "mean_tls": mean_tls,
                "status": test_status
            })
    
    core_tls_df = pd.DataFrame(results_core_tls)
    core_tls_df.set_index(["solidity_type", "signature"], inplace=True)
    
    # ===== Part 2: Pairwise solidity comparisons within each architecture =====
    results_pairwise = []
    
    for arch in ["Core", "TLS"]:
        df_arch = df[df["architecture_type"] == arch].copy()
        
        # Generate all pairwise comparisons
        from itertools import combinations
        for sol1, sol2 in combinations(solidity_order, 2):
            df_sol1 = df_arch[df_arch["solidity_type"] == sol1]
            df_sol2 = df_arch[df_arch["solidity_type"] == sol2]
            
            for sig in sig_cols:
                x1 = df_sol1[sig].dropna()
                x2 = df_sol2[sig].dropna()
                
                n1 = len(x1)
                n2 = len(x2)
                
                if n1 >= min_samples and n2 >= min_samples:
                    try:
                        _, pval = mannwhitneyu(x1, x2, alternative="two-sided")
                        test_status = "ok"
                    except Exception as e:
                        pval = np.nan
                        test_status = f"error: {str(e)}"
                else:
                    pval = np.nan
                    test_status = "insufficient"
                
                results_pairwise.append({
                    "architecture_type": arch,
                    "signature": sig,
                    "comparison": f"{sol1}_vs_{sol2}",
                    "pval": pval,
                    "n1": n1,
                    "n2": n2,
                    "status": test_status
                })
    
    pairwise_df = pd.DataFrame(results_pairwise)
    pairwise_df.set_index(["architecture_type", "signature", "comparison"], inplace=True)
    
    return core_tls_df, pairwise_df


def plot_signature_boxplots_by_solidity(
    df, 
    core_tls_stats_df,
    pairwise_stats_df,
    sig_prefix="sig:",
    sig_suffix="_z",
    solidity_order=None, 
    figsize=(12, 6),
    pval_threshold=0.05,
    show_ns=True,
    output_dir="figures/stats"
):
    """
    Plots boxplots comparing Core vs TLS across solidity types for each signature.
    Creates one plot per signature with solidity types on x-axis.
    Shows both Core vs TLS comparisons and pairwise solidity comparisons.
    
    Parameters
    ----------
    df : pd.DataFrame
        Data containing signature columns and grouping variables
    core_tls_stats_df : pd.DataFrame
        Statistics for Core vs TLS comparisons
    pairwise_stats_df : pd.DataFrame
        Statistics for pairwise solidity comparisons
    sig_prefix : str
        Prefix for signature columns
    solidity_order : list or None
        Order of solidity groups
    figsize : tuple
        Figure size (width, height)
    pval_threshold : float
        P-value threshold for significance
    show_ns : bool
        Whether to show "n.s." for non-significant comparisons
    output_dir : str
        Directory to save figures
    """
    import os
    
    # Get signature columns
    sig_cols = [c for c in df.columns if c.startswith(sig_prefix) and c.endswith(sig_suffix)]
    
    if len(sig_cols) == 0:
        print(f"Warning: No signature columns found starting with '{sig_prefix}'")
        return
    
    if solidity_order is None:
        solidity_order = sorted(df["solidity_type"].dropna().unique())
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Filter data once
    df_filtered = df[df["architecture_type"].isin(["Core", "TLS"])].copy()
    
    # Create one figure per signature
    for sig in sig_cols:
        print(f"Plotting {sig}...")
        
        fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)
        
        # Left panel: Core, Right panel: TLS
        for ax_idx, arch_type in enumerate(["Core", "TLS"]):
            ax = axes[ax_idx]
            
            # Prepare data for this architecture type across all solidity types
            data_to_plot = []
            labels = []
            sample_sizes = []
            
            for sol in solidity_order:
                df_subset = df_filtered[
                    (df_filtered["solidity_type"] == sol) & 
                    (df_filtered["architecture_type"] == arch_type)
                ][sig].dropna()
                
                data_to_plot.append(df_subset.values)
                labels.append(sol)
                sample_sizes.append(len(df_subset))
            
            # Check if we have any data
            if all(len(d) == 0 for d in data_to_plot):
                ax.text(0.5, 0.5, f"No {arch_type} data", 
                       ha="center", va="center", transform=ax.transAxes)
                ax.set_title(arch_type, fontsize=12, fontweight='bold')
                ax.set_axis_off()
                continue
            
            # Plot boxplots
            try:
                positions = np.arange(1, len(labels) + 1)
                bp = ax.boxplot(
                    data_to_plot,
                    labels=labels,
                    positions=positions,
                    patch_artist=True,
                    widths=0.6
                )
                
                # Color the boxes
                colors = ['#e0e0e0', '#a6cee3', '#1f78b4']  # Light to dark
                for patch, color in zip(bp['boxes'], colors[:len(labels)]):
                    patch.set_facecolor(color)
                    patch.set_alpha(0.7)
                
                # Add stripplot for individual points
                for i, data in enumerate(data_to_plot, start=1):
                    if len(data) > 0:
                        y = data
                        x = np.random.normal(i, 0.04, size=len(y))
                        ax.scatter(x, y, alpha=0.4, s=15, color='black')
                
                # Get y-axis limits for positioning annotations
                y_min, y_max = ax.get_ylim()
                y_range = y_max - y_min
                
                # Add sample sizes at the bottom of each box
                for i, (label, n) in enumerate(zip(labels, sample_sizes), start=1):
                    ax.text(i, y_min - 0.05 * y_range, f"n={n}", 
                           ha="center", va="top", fontsize=8, 
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                                    edgecolor='gray', alpha=0.8))
                
                # Add pairwise comparison p-values ABOVE the boxes
                # Comparisons: 1 vs 2, 2 vs 3, 1 vs 3
                comparisons = [
                    (0, 1, "Normal_vs_Non-Solid", 0.15),
                    (1, 2, "Non-Solid_vs_Solid", 0.10),
                    (0, 2, "Normal_vs_Solid", 0.20)
                ]
                
                for idx1, idx2, comp_name, y_offset in comparisons:
                    if idx1 >= len(labels) or idx2 >= len(labels):
                        continue
                    
                    # Look up p-value
                    pval = np.nan
                    if (arch_type, sig, comp_name) in pairwise_stats_df.index:
                        row = pairwise_stats_df.loc[(arch_type, sig, comp_name)]
                        pval = row["pval"]
                        status = row["status"]
                    
                    # Format p-value label
                    if pd.notna(pval) and status == "ok":
                        if pval < 0.001:
                            label = "***"
                            color = "red"
                        elif pval < 0.01:
                            label = "**"
                            color = "red"
                        elif pval < pval_threshold:
                            label = "*"
                            color = "red"
                        else:
                            if show_ns:
                                label = "n.s."
                                color = "gray"
                            else:
                                label = ""
                                color = "gray"
                    else:
                        if show_ns:
                            label = "n.s."
                            color = "gray"
                        else:
                            label = ""
                            color = "gray"
                    
                    # Draw significance bar
                    if label:
                        x1 = positions[idx1]
                        x2 = positions[idx2]
                        y_pos = y_max + y_offset * y_range
                        
                        ax.plot([x1, x1, x2, x2], 
                               [y_pos - 0.02*y_range, y_pos, y_pos, y_pos - 0.02*y_range], 
                               color='black', linewidth=1)
                        ax.text((x1 + x2) / 2, y_pos + 0.01*y_range, label, 
                               ha="center", va="bottom", fontsize=9, color=color,
                               fontweight='bold' if color == 'red' else 'normal')
                
            except Exception as e:
                ax.text(0.5, 0.5, f"Plot error:\n{str(e)}", 
                       ha="center", va="center", transform=ax.transAxes, fontsize=8)
            
            # Formatting
            ax.set_title(arch_type, fontsize=14, fontweight='bold', pad=20)
            ax.set_xlabel("Solidity Type", fontsize=11)
            if ax_idx == 0:
                ax.set_ylabel("Score", fontsize=11)
            ax.tick_params(labelsize=10)
            ax.grid(axis='y', alpha=0.3, linestyle='--')
            
            # Extend y-axis to accommodate annotations
            y_min, y_max = ax.get_ylim()
            y_range = y_max - y_min
            ax.set_ylim(y_min - 0.1 * y_range, y_max + 0.35 * y_range)
        
        # Add Core vs TLS comparison p-values as text in the figure
        pval_text = "Core vs TLS:  "
        for sol in solidity_order:
            if (sol, sig) in core_tls_stats_df.index:
                row = core_tls_stats_df.loc[(sol, sig)]
                pval = row["pval"]
                status = row["status"]
                
                if pd.notna(pval) and status == "ok":
                    if pval < 0.001:
                        pval_text += f"{sol}: p<0.001  "
                    elif pval < pval_threshold:
                        pval_text += f"{sol}: p={pval:.3f}  "
                    else:
                        if show_ns:
                            pval_text += f"{sol}: n.s.  "
                else:
                    if show_ns:
                        pval_text += f"{sol}: n.s.  "
        
        fig.text(0.5, 0.02, pval_text, ha="center", va="bottom", fontsize=10,
                bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow', alpha=0.8))
        
        # Main title
        sig_clean = sig.replace(sig_prefix, "")
        plt.suptitle(f"{sig_clean}", fontsize=14, fontweight='bold', y=0.98)
        
        plt.tight_layout(rect=[0, 0.05, 1, 0.96])
        
        # Save figure with cleaned filename
        filename = sig.replace(sig_prefix, "").replace("/", "_").replace(":", "_")
        filepath = os.path.join(output_dir, f"{filename}.png")
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved to {filepath}")


# Run statistics (now returns TWO dataframes)
stats_core_tls, stats_pairwise = compare_core_tls_by_solidity(
    df, 
    sig_prefix="sig:",
    sig_suffix="_z",
    solidity_order=["Normal", "Non-Solid", "Solid"],
    min_samples=3
)

# Plot (now takes both dataframes)
plot_signature_boxplots_by_solidity(
    df,
    stats_core_tls,
    stats_pairwise,
    sig_prefix="sig:",
    sig_suffix="_z",
    solidity_order=["Normal", "Non-Solid", "Solid"],
    output_dir="figures/stats"
)







# P1_tumor_samples = ['PT01-2_TUM', 'PT01-4_TUM', 'PT01-6_TUM']
# p1 = adata[adata.obs['library_id'].isin(P1_tumor_samples)].copy()
p1 = adata.copy()

# choose the columns you actually have (z-scored recommended if you centered them)
tumor_col = 'sig:CRC_Tumor/Epithelial_Core'         # or 'sig:CRC_Tumor/Epithelial_Core_z'
cycle_col = 'sig:Proliferation/Cell_Cycle_z'        # ensure this column exists

# simple thresholds (adjust as needed)
p1.obs['tumor'] = p1.obs[tumor_col] > 0

p1.obs['proliferative_tumor'] = (
    (p1.obs[cycle_col] > 0) &
    (p1.obs[tumor_col] > 0)
)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- Preconditions ---
assert "library_id" in p1.obs.columns, "p1.obs['library_id'] is required"
for col in ["tumor", "proliferative_tumor"]:
    assert col in p1.obs.columns, f"Missing boolean column: {col}"

# --- Aggregate per library_id ---
grp = p1.obs.groupby("library_id", observed=True)

summary = pd.DataFrame({
    "n_spots": grp.size(),
    "n_tumor": grp["tumor"].sum(),
    "n_prolif_tumor": grp["proliferative_tumor"].sum(),
})

summary["prop_tumor_over_all"] = summary["n_tumor"] / summary["n_spots"]
summary["prop_prolif_over_tumor"] = np.where(
    summary["n_tumor"] > 0,
    summary["n_prolif_tumor"] / summary["n_tumor"],
    0.0,
)

# Order by your P1 samples
order = [lib for lib in p1.obs["library_id"].unique().tolist()]
summary = summary.reindex(order)

# --- Plot separate barplots ---
fig, axes = plt.subplots(2, 1, figsize=(3, 4))

def annotate_bars(ax, values):
    for i, v in enumerate(values):
        ax.annotate(f"{v*100:.2f}%", (i, v),
                    xytext=(0, 3), textcoords="offset points",
                    ha="center", va="bottom", fontsize=9, clip_on=False)

# 1. Tumor / All
axes[0].bar(summary.index, summary["prop_tumor_over_all"], color="steelblue")
axes[0].set_ylabel("Tumor / All")
axes[0].set_ylim(0, 1.15)
axes[0].set_title("")
axes[0].grid(False)
annotate_bars(axes[0], summary["prop_tumor_over_all"])
axes[0].set_xticks([0, 1,2], ['BL', 'PRT', 'Surgery'])

# 2. Proliferative Tumor / Tumor
axes[1].bar(summary.index, summary["prop_prolif_over_tumor"], color="steelblue")
axes[1].set_ylabel("Prolif. Tumor / Tumor")
axes[1].set_ylim(0, 1.15)
axes[1].grid(False)
axes[1].set_title("")
annotate_bars(axes[1], summary["prop_prolif_over_tumor"])

axes[1].set_xlabel("")
axes[1].set_xticks([0, 1,2],['BL', 'PRT', 'Surgery'])
sns.despine()
plt.tight_layout()
plt.show()








import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.spatial import cKDTree

# --- Preconditions ---
assert "library_id" in p1.obs.columns
assert "tumor" in p1.obs.columns              # boolean mask per spot
assert "sig:Myeloid/Macrophage_Repair_z" in p1.obs.columns
assert "spatial" in p1.obsm, "p1.obsm['spatial'] missing"

# define epithelial repair positivity if not already present
if "macrophage_repair" not in p1.obs.columns:
    p1.obs["macrophage_repair"] = (p1.obs["sig:Myeloid/Macrophage_Repair_z"] > 0)

# ----------- Helper: periphery by pixel distance per library -----------
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

def define_boundary_band(
    adata_sub,
    periphery_px: float = 80.0,
    obsm_key: str = "spatial",
    tumor_col: str = "tumor",
    mode: str = "both",  # "non_tumor_only", "tumor_only", or "both"
) -> pd.Series:
    """
    Return a boolean Series indexed like adata_sub.obs where True marks spots
    near the tumor boundary within `periphery_px` pixels.

    Logic
      - Build a KDTree on tumor spots and query all non-tumor distances.
      - Build a KDTree on non-tumor spots and query all tumor distances.
      - Modes:
          "non_tumor_only": non-tumor spots with dist_to_tumor <= periphery_px
          "tumor_only"    : tumor spots with dist_to_non_tumor <= periphery_px
          "both"          : union of the two sets.
    Computed per library_id to avoid cross-slide mixing.
    """
    assert obsm_key in adata_sub.obsm, f"{obsm_key} not found in adata_sub.obsm"
    assert tumor_col in adata_sub.obs.columns, f"{tumor_col} not in adata_sub.obs"
    assert "library_id" in adata_sub.obs.columns, "library_id missing in obs"
    assert mode in {"non_tumor_only", "tumor_only", "both"}, "invalid mode"

    XY_all = adata_sub.obsm[obsm_key]
    is_tumor_all = adata_sub.obs[tumor_col].astype(bool).values
    libs = adata_sub.obs["library_id"].values

    band = np.zeros(adata_sub.n_obs, dtype=bool)

    for lib in np.unique(libs):
        sel = np.where(libs == lib)[0]
        if sel.size == 0:
            continue

        XY = XY_all[sel]
        is_tumor = is_tumor_all[sel]

        # Split by class
        idx_tumor = np.where(is_tumor)[0]
        idx_ntum  = np.where(~is_tumor)[0]

        # If either class is absent, no boundary can be defined for that library
        if idx_tumor.size == 0 or idx_ntum.size == 0:
            continue

        # Distances: non-tumor -> nearest tumor
        tree_tumor = cKDTree(XY[idx_tumor])
        d_ntum_to_tum, _ = tree_tumor.query(XY[idx_ntum], k=1)

        # Distances: tumor -> nearest non-tumor
        tree_ntum = cKDTree(XY[idx_ntum])
        d_tum_to_ntum, _ = tree_ntum.query(XY[idx_tumor], k=1)

        mark = np.zeros(sel.size, dtype=bool)

        if mode in {"non_tumor_only", "both"}:
            mark[idx_ntum[d_ntum_to_tum <= periphery_px]] = True
        if mode in {"tumor_only", "both"}:
            mark[idx_tumor[d_tum_to_ntum <= periphery_px]] = True

        band[sel] = mark

    return pd.Series(band, index=adata_sub.obs.index, name=f"boundary_{mode}")

# # choose your pixel distance threshold
periphery_px = 10
p1.obs["periphery"] = define_boundary_band(p1, periphery_px=2, mode="both")

# --- Aggregate per library_id ---
grp = p1.obs.groupby("library_id", observed=True)

# counts
summary = pd.DataFrame({
    "n_spots": grp.size(),
    "n_tumor": grp["tumor"].sum(),
    "n_repair_all": grp["macrophage_repair"].sum(),
    "n_periphery": grp["periphery"].sum(),
    "n_repair_in_tumor": (grp.apply(lambda df: (df["macrophage_repair"] & df["tumor"]).sum())),
    "n_repair_in_periphery": (grp.apply(lambda df: (df["macrophage_repair"] & df["periphery"]).sum())),
})

# proportions
summary["prop_repair_over_tumor"] = np.where(
    summary["n_tumor"] > 0,
    summary["n_repair_in_tumor"] / summary["n_tumor"],
    0.0,
)
summary["prop_repair_over_periphery"] = np.where(
    summary["n_periphery"] > 0,
    summary["n_repair_in_periphery"] / summary["n_periphery"],
    0.0,
)

# keep your order
order = p1.obs["library_id"].unique().tolist()
summary = summary.reindex(order)

# optional short x tick labels
xtick_labels = ["BL", "PRT", "Surgery"] if len(summary.index) == 3 else summary.index.tolist()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# --- Plot: two separate barplots, independent axes ---
fig, axes = plt.subplots(2, 1, figsize=(3, 4), sharex=False)
def annotate_bars(ax, values):
    for i, v in enumerate(values):
        ax.annotate(f"{v*100:.2f}%", (i, v),
                    xytext=(0, 3), textcoords="offset points",
                    ha="center", va="bottom", fontsize=9, clip_on=False)

# 1) Repair / Tumor
vals1 = summary["prop_repair_over_tumor"].values
axes[0].bar(range(len(summary.index)), vals1, color="steelblue")
axes[0].set_ylim(0, 1.15)  # ymin=0, ymax automatic
axes[0].set_ylabel("Repair / Tumor")
axes[0].grid(False)
axes[0].set_xticks(range(len(summary.index)))
axes[0].set_xticklabels(xtick_labels)
annotate_bars(axes[0], vals1)

# 2) Repair / Periphery
vals2 = summary["prop_repair_over_periphery"].values
axes[1].bar(range(len(summary.index)), vals2, color="steelblue")
axes[1].set_ylim(0, None)  # ymin=0, ymax automatic
axes[1].set_ylabel("Repair / Periphery")
axes[1].grid(False)
axes[1].set_xticks(range(len(summary.index)))
axes[1].set_xticklabels(xtick_labels)
annotate_bars(axes[1], vals2)

sns.despine()
plt.tight_layout()
plt.show()


fig, axes = plt.subplots(1,3)
for i, library_id in enumerate(p1.obs['library_id'].unique()):
    p1_ = p1[p1.obs['library_id'] == library_id].copy()
    sc.pl.spatial(p1_, color = 'periphery', library_id = library_id, frameon = False, show = False, ax = axes[i], title = '', legend_loc = None, palette = ['lightgray', 'orange'])
plt.show()





import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc

# --- thresholds (edit if you prefer different cutoffs) ---
PROLIF_THR = 0.0
REPAIR_THR = 0.0

def _to_categorical_zone(df, ref_mask, cmp_mask):
    """
    Return a categorical pd.Series with values in {"rest","reference","compare"}.
    'compare' takes precedence over 'reference' if both True.
    """
    out = np.full(df.shape[0], "rest", dtype=object)
    out[ref_mask] = "reference"
    out[cmp_mask] = "compare"    # overwrite if overlap
    return out

def _attach_colors(adata_sub, col):
    """
    Bind explicit colors to a categorical obs column so that:
      rest -> transparent, reference -> yellow, compare -> purple
    """
    cats = list(adata_sub.obs[col].astype("category").cat.categories)
    lut = {"rest": "#lightblue", "reference": "yellow", "compare": "purple"}
    adata_sub.uns[f"{col}_colors"] = [lut[c] for c in cats]

def _plot_zone(ax, adata_sub, col, library_id):
    _attach_colors(adata_sub, col)
    sc.pl.spatial(
        adata_sub,
        color=col,
        library_id=library_id,
        frameon=False,
        show=False,
        ax=ax,
        title="",
        legend_loc=None,   # no legend
    )

# ---------- build zone columns once ----------
# Tumor vs Rest
p1.obs["zone_tumor"] = _to_categorical_zone(
    p1.obs,
    ref_mask=p1.obs["tumor"].values,
    cmp_mask=np.zeros(p1.n_obs, dtype=bool)
)

# Prolif Tumor vs Tumor vs Rest
prolif_tumor = (p1.obs["sig:Proliferation/Cell_Cycle_z"].values > PROLIF_THR) & p1.obs["tumor"].values
p1.obs["zone_prolifTumor"] = _to_categorical_zone(
    p1.obs,
    ref_mask=p1.obs["tumor"].values,
    cmp_mask=prolif_tumor
)

# Repair Tumor vs Tumor vs Rest
repair_pos = (p1.obs["sig:Stromal/Epithelial_Repair_z"].values > REPAIR_THR)
repair_tumor = repair_pos & p1.obs["tumor"].values
p1.obs["zone_repairTumor"] = _to_categorical_zone(
    p1.obs,
    ref_mask=p1.obs["tumor"].values,
    cmp_mask=repair_tumor
)

# Repair Periphery vs Periphery vs Rest
# expects p1.obs["periphery_both"] already computed
repair_periphery = repair_pos & p1.obs["periphery"].values
p1.obs["zone_repairPeriphery"] = _to_categorical_zone(
    p1.obs,
    ref_mask=p1.obs["periphery"].values,
    cmp_mask=repair_periphery
)

# ---------- plot loop ----------
for library_id in p1.obs["library_id"].unique():
    p1_ = p1[p1.obs["library_id"] == library_id].copy()
    fig, axes = plt.subplots(1, 4, figsize=(12, 4))
    cols = ["zone_tumor", "zone_prolifTumor", "zone_repairTumor", "zone_repairPeriphery"]

    for ax, col in zip(axes, cols):
        _plot_zone(ax, p1_, col, library_id)

    # optional: small spacing and file naming
    plt.subplots_adjust(wspace=0.02)
    plt.suptitle(library_id, y=0.98)
    plt.savefig(f"figures/spatial/{library_id}_zones.png", dpi=300, bbox_inches="tight")
    plt.close()

zone_cols = {
    "Tumor vs All": "zone_tumor",
    "Proliferative Tumor vs Tumor": "zone_prolifTumor",
    "Repair Tumor vs Tumor": "zone_repairTumor",
    "Repair Periphery vs Periphery": "zone_repairPeriphery",
}

# --- plot loop per zone ---
for zone_name, col in zone_cols.items():
    libs = p1.obs["library_id"].unique().tolist()
    fig, axes = plt.subplots(1, len(libs), figsize=(4*len(libs), 4))

    if len(libs) == 1:
        axes = [axes]

    for ax, library_id in zip(axes, libs):
        p1_ = p1[p1.obs["library_id"] == library_id].copy()
        _plot_zone(ax, p1_, col, library_id)
        ax.set_title('', fontsize=9)

    plt.subplots_adjust(wspace=0.02)
    plt.suptitle(zone_name, y=0.95, fontsize=12)
    plt.savefig(f"figures/spatial/{zone_name}_allLibraries.png", dpi=300, bbox_inches="tight")
    plt.close()




import numpy as np
import pandas as pd
import scanpy as sc

# --- thresholds ---
PROLIF_THR = 0.0
REPAIR_THR = 0.0

def to_zone_labels(n, ref_mask, cmp_mask):
    """
    Build categorical labels in {"rest","reference","compare"} where
    'reference' overwrites 'compare' in overlapping spots. This ensures
    the subset stays orange over the broader background in gray.
    """
    out = np.full(n, "rest", dtype=object)
    out[cmp_mask] = "compare"
    out[ref_mask] = "reference"  # reference wins
    return pd.Categorical(out, categories=["rest", "reference", "compare"])

def attach_zone_colors(adata_sub, col):
    """
    Bind colors per category: rest=transparent, reference=orange, compare=lightgray.
    """
    cats = list(adata_sub.obs[col].cat.categories)
    lut = {"rest": "#00000000", "reference": "orange", "compare": "lightgray"}
    adata_sub.uns[f"{col}_colors"] = [lut[c] for c in cats]

# sanity
assert "tumor" in p1.obs.columns
assert "sig:Proliferation/Cell_Cycle_z" in p1.obs.columns
assert "sig:Stromal/Epithelial_Repair_z" in p1.obs.columns
assert "periphery" in p1.obs.columns  # from boundary band helper

# convenience masks
tumor = p1.obs["tumor"].values.astype(bool)
prolif_pos = (p1.obs["sig:Proliferation/Cell_Cycle_z"].values > PROLIF_THR)
repair_pos  = (p1.obs["sig:Stromal/Epithelial_Repair_z"].values > REPAIR_THR)
periphery = p1.obs["periphery"].values.astype(bool)

# 1) Tumor vs Rest  (only reference shown; rest is transparent)
p1.obs["zone_tumor"] = to_zone_labels(
    p1.n_obs,
    ref_mask=tumor,
    cmp_mask=np.zeros(p1.n_obs, dtype=bool),
)

# 2) Proliferative Tumor vs Tumor vs Rest
# reference = proliferative tumor (orange), compare = tumor (lightgray)
prolif_tumor = prolif_pos & tumor
p1.obs["zone_prolifTumor"] = to_zone_labels(
    p1.n_obs,
    ref_mask=prolif_tumor,
    cmp_mask=tumor,
)

# 3) Repair Tumor vs Tumor vs Rest
repair_tumor = repair_pos & tumor
p1.obs["zone_repairTumor"] = to_zone_labels(
    p1.n_obs,
    ref_mask=repair_tumor,
    cmp_mask=tumor,
)

# # 4) Repair Periphery vs Periphery vs Rest
# repair_periphery = repair_pos & periphery
# p1.obs["zone_repairPeriphery"] = to_zone_labels(
#     p1.n_obs,
#     ref_mask=repair_periphery,
#     cmp_mask=periphery,
# )

# # attach colors for each zone column just once
# for col in ["zone_tumor", "zone_prolifTumor", "zone_repairTumor", "zone_repairPeriphery"]:
#     p1.obs[col] = p1.obs[col].astype("category")
#     attach_zone_colors(p1, col)
