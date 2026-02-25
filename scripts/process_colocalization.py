"""
Spot Process Colocalization Analysis: Measure spatial co-occurrence patterns of processes.

This script analyzes spatial co-occurrence patterns of gene signature scores (processes)
across spots using multiple modular methods. Each analysis can be enabled/disabled via
the ANALYSIS_CONFIG dictionary.

Available analyses:
1. Pearson correlation of processes across spots (using DataFrame.corr())
2. Moran's I for spatial autocorrelation of processes (using sq.gr.spatial_autocorr)
3. Thresholded neighborhood enrichment using squidpy.gr.nhood_enrichment
   after thresholding processes by low/high (-1, +1)

Output: Heatmaps showing colocalization patterns and identify process pairs that are
colocalized (positive association) or anti-colocalized (mutual exclusion).
"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import squidpy as sq
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

sc.settings.set_figure_params(dpi=300, dpi_save=400)

# ============================================================================
# Configuration: Enable/Disable Analyses
# ============================================================================

ANALYSIS_CONFIG = {
    # Set to True to run each analysis, False to skip
    "pearson_correlation": True,  # Pearson correlation matrix across spots
    "morans_i": False,  # Moran's I spatial autocorrelation
    "neighborhood_enrichment": False,  # Thresholded neighborhood enrichment (slower)
    "volcano_plots": True,  # Volcano plots for Normal vs Non-Solid vs Solid
}

# Analysis parameters
THRESHOLD_LOW = -1.0  # Threshold for low classification (z-scores)
THRESHOLD_HIGH = 1.0  # Threshold for high classification (z-scores)
N_PERMS = 1000  # Number of permutations for neighborhood enrichment

# Volcano plot signature selection
# Set to None to use all signatures, or provide lists to include/exclude
# Supports both exact signature names (with or without 'sig:' prefix and '_z' suffix)
# and pattern matching (substring search, case-insensitive)
VOLCANO_SIGNATURES_INCLUDE = None  # List of signatures to include (None = all)
VOLCANO_SIGNATURES_EXCLUDE = [
    "sig:Immune_Myeloid/Classical_Monocyte",
    "sig:Immune_Myeloid/NonClassical_Monocyte",
    "sig:Immune_Myeloid/Alveolar_Macrophage",
    "sig:Immune_Myeloid/MonoDerived_Macrophage",
    "sig:Immune_Myeloid/M1_Macrophage",
    "sig:Immune_Myeloid/M2_Macrophage",
]  # List of signatures to exclude (None = none)

# Examples:
# VOLCANO_SIGNATURES_INCLUDE = ['sig:Tumor_Cells-EMT_Tumor_z', 'Hypoxia']  # Exact or pattern
# VOLCANO_SIGNATURES_EXCLUDE = ['TLS', 'Immune_Lymphoid']  # Exclude patterns

# ============================================================================
# Helper Functions
# ============================================================================


def get_signature_columns(adata, prefix="sig:"):
    """Get all signature columns from adata.obs."""
    sig_cols = [col for col in adata.obs.columns if col.startswith(prefix) and col.endswith("_z")]
    return sorted(sig_cols)


def filter_signatures(signature_list, include=None, exclude=None):
    """
    Filter signatures based on include/exclude criteria.

    Parameters
    ----------
    signature_list : list
        List of signature column names (e.g., 'sig:Tumor_Cells-EMT_Tumor_z')
    include : list or None
        List of signatures to include. Can be:
        - Exact matches (with or without 'sig:' prefix and '_z' suffix)
        - Patterns (substring search, case-insensitive)
        - None to include all
    exclude : list or None
        List of signatures to exclude. Same format as include.
        - None to exclude none

    Returns
    -------
    list
        Filtered list of signatures
    """
    if include is None and exclude is None:
        return signature_list

    filtered = []

    # Normalize signature names for matching (remove prefix/suffix)
    def normalize_sig(sig):
        """Remove prefix and suffix for comparison."""
        normalized = sig
        if normalized.startswith("sig:"):
            normalized = normalized[4:]
        if normalized.endswith("_z"):
            normalized = normalized[:-2]
        return normalized.lower()

    # Normalize all signatures
    sig_normalized = {sig: normalize_sig(sig) for sig in signature_list}

    for sig in signature_list:
        sig_norm = sig_normalized[sig]
        include_match = False
        exclude_match = False

        # Check include list
        if include is not None:
            for pattern in include:
                # Normalize pattern
                if pattern.startswith("sig:"):
                    pattern_norm = normalize_sig(pattern)
                else:
                    pattern_norm = pattern.lower()

                # Check exact match (normalized) or substring match (pattern in signature)
                if pattern_norm == sig_norm or pattern_norm in sig_norm:
                    include_match = True
                    break
        else:
            # No include list means include all
            include_match = True

        # Check exclude list
        if exclude is not None:
            for pattern in exclude:
                # Normalize pattern
                if pattern.startswith("sig:"):
                    pattern_norm = normalize_sig(pattern)
                else:
                    pattern_norm = pattern.lower()

                # Check exact match (normalized) or substring match (pattern in signature)
                if pattern_norm == sig_norm or pattern_norm in sig_norm:
                    exclude_match = True
                    break

        # Include if matches include criteria and doesn't match exclude
        if include_match and not exclude_match:
            filtered.append(sig)

    return filtered


def threshold_signature(adata, sig_col, threshold_low=-1.0, threshold_high=1.0):
    """
    Threshold signature scores into low, medium, high categories.

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    sig_col : str
        Signature column name
    threshold_low : float
        Threshold for low category (default: -1.0)
    threshold_high : float
        Threshold for high category (default: 1.0)

    Returns
    -------
    np.ndarray
        Categorical array: 'low', 'medium', 'high'
    """
    scores = adata.obs[sig_col].values
    categories = np.full(len(scores), "medium", dtype=object)
    categories[scores <= threshold_low] = "low"
    categories[scores >= threshold_high] = "high"
    return categories


# ============================================================================
# Main Analysis Functions
# ============================================================================

# ============================================================================
# Modular Analysis Functions
# ============================================================================


def analyze_pearson_correlation(adata, sig_columns, output_dir):
    """
    Analyze Pearson correlation between processes across spots.

    Parameters
    ----------
    adata : AnnData
        Annotated data object with signature scores
    sig_columns : list
        List of signature columns to analyze
    output_dir : str
        Output directory for results

    Returns
    -------
    pd.DataFrame
        Correlation matrix
    """
    print("\n" + "=" * 70)
    print("ANALYSIS 1: PEARSON CORRELATION")
    print("=" * 70)

    # Extract signature scores as DataFrame
    print("\nExtracting signature scores...")
    sig_df = adata.obs[sig_columns].copy()

    # Remove signatures with too many NaN values
    valid_sigs = sig_df.columns[sig_df.isna().sum() < sig_df.shape[0] * 0.5].tolist()
    sig_df = sig_df[valid_sigs]
    print(f"   Using {len(valid_sigs)} signatures with <50% missing values")

    # Compute correlation matrix
    print("\nComputing Pearson correlation matrix...")
    corr_matrix = sig_df.corr(method="pearson")

    # Generate standard heatmap
    print("\nGenerating correlation heatmap...")
    fig, ax = plt.subplots(figsize=(max(12, len(valid_sigs) * 0.5), max(10, len(valid_sigs) * 0.5)))
    sns.heatmap(
        corr_matrix,
        cmap="RdBu_r",
        center=0,
        vmin=-1,
        vmax=1,
        square=True,
        cbar_kws={"label": "Pearson Correlation"},
        ax=ax,
        fmt=".2f",
        annot=False,
    )
    ax.set_title(
        "Process Colocalization: Pearson Correlation Matrix", fontsize=14, fontweight="bold"
    )
    plt.tight_layout()
    plt.savefig(
        "figures/process_colocalization/correlation_heatmap.pdf", bbox_inches="tight", dpi=300
    )
    plt.savefig(
        "figures/process_colocalization/correlation_heatmap.png", bbox_inches="tight", dpi=300
    )
    plt.close()
    print("   ✅ Correlation heatmap saved")

    # Generate clustermap (hierarchically clustered heatmap)
    print("\nGenerating correlation clustermap...")
    g = sns.clustermap(
        corr_matrix,
        cmap="RdBu_r",
        center=0,
        vmin=-1,
        vmax=1,
        figsize=(max(14, len(valid_sigs) * 0.6), max(12, len(valid_sigs) * 0.6)),
        linewidths=0.5,
        cbar_kws={"label": "Pearson Correlation"},
        dendrogram_ratio=0.1,
        cbar_pos=(0.02, 0.83, 0.03, 0.12),
        row_cluster=True,
        col_cluster=True,
        yticklabels=True,
        xticklabels=True,
        method="average",  # Linkage method
        metric="euclidean",  # Distance metric
    )

    # Set title
    g.fig.suptitle(
        "Process Colocalization: Pearson Correlation (Hierarchically Clustered)",
        fontsize=14,
        fontweight="bold",
        y=0.995,
    )

    # Adjust labels
    g.ax_heatmap.set_xlabel("Process Signatures", fontsize=11)
    g.ax_heatmap.set_ylabel("Process Signatures", fontsize=11)

    plt.savefig(
        "figures/process_colocalization/correlation_clustermap.pdf", bbox_inches="tight", dpi=300
    )
    plt.savefig(
        "figures/process_colocalization/correlation_clustermap.png", bbox_inches="tight", dpi=300
    )
    plt.close()
    print("   ✅ Correlation clustermap saved")

    # Save results
    corr_matrix.to_csv(os.path.join(output_dir, "correlation_matrix.csv"))
    print(f"   ✅ Results saved to: {output_dir}/correlation_matrix.csv")

    # Print summary
    corr_pairs = []
    for i, sig1 in enumerate(valid_sigs):
        for sig2 in valid_sigs[i + 1 :]:
            r = corr_matrix.loc[sig1, sig2]
            if not np.isnan(r):
                corr_pairs.append((sig1, sig2, r))

    corr_pairs_sorted = sorted(corr_pairs, key=lambda x: x[2], reverse=True)

    print("\nTop 5 colocalized pairs (highest correlation):")
    for sig1, sig2, r in corr_pairs_sorted[:5]:
        print(f"   {sig1} <-> {sig2}: r = {r:.3f}")

    print("\nTop 5 anti-colocalized pairs (lowest correlation):")
    for sig1, sig2, r in corr_pairs_sorted[-5:]:
        print(f"   {sig1} <-> {sig2}: r = {r:.3f}")

    return corr_matrix


def analyze_morans_i(adata, sig_columns, output_dir):
    """
    Analyze Moran's I spatial autocorrelation for each process.

    Uses manual computation via spatial neighbors graph since spatial_autocorr
    expects genes in var_names, not signature scores in obs.

    Parameters
    ----------
    adata : AnnData
        Annotated data object with signature scores
    sig_columns : list
        List of signature columns to analyze
    output_dir : str
        Output directory for results

    Returns
    -------
    pd.DataFrame
        DataFrame with Moran's I results
    """
    print("\n" + "=" * 70)
    print("ANALYSIS 2: MORAN'S I SPATIAL AUTOCORRELATION")
    print("=" * 70)

    # Extract signature scores
    sig_df = adata.obs[sig_columns].copy()
    valid_sigs = sig_df.columns[sig_df.isna().sum() < sig_df.shape[0] * 0.5].tolist()
    print(f"\nAnalyzing {len(valid_sigs)} signatures")

    # Compute spatial neighbors
    print("\nComputing spatial neighbors...")
    if "spatial_neighbors" not in adata.obsp:
        sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)

    # Compute Moran's I manually using spatial neighbors graph
    print("\nComputing Moran's I for each signature...")
    moran_results = {}

    # Get connectivity matrix
    if "spatial_connectivities" in adata.obsp:
        W = adata.obsp["spatial_connectivities"]
        if hasattr(W, "toarray"):
            W_dense = W.toarray()
        else:
            W_dense = W
    else:
        print("   Error: spatial_connectivities not found")
        return pd.DataFrame()

    n = adata.n_obs

    for sig_col in tqdm(valid_sigs, desc="   Computing Moran's I"):
        scores = adata.obs[sig_col].values

        # Remove NaN values
        valid_mask = ~np.isnan(scores)
        if valid_mask.sum() < 3:
            moran_results[sig_col] = {"I": np.nan, "pval": np.nan}
            continue

        scores_clean = scores[valid_mask]
        W_sub = W_dense[valid_mask, :][:, valid_mask]

        n_valid = len(scores_clean)
        if n_valid < 3:
            moran_results[sig_col] = {"I": np.nan, "pval": np.nan}
            continue

        # Compute mean
        x_mean = np.mean(scores_clean)

        # Compute deviations
        x_centered = scores_clean - x_mean

        # Compute sum of squared deviations
        sum_sq_dev = np.sum(x_centered**2)
        if sum_sq_dev == 0:
            moran_results[sig_col] = {"I": np.nan, "pval": np.nan}
            continue

        # Compute W (sum of all weights)
        W_sum = np.sum(W_sub)
        if W_sum == 0:
            moran_results[sig_col] = {"I": np.nan, "pval": np.nan}
            continue

        # Compute Moran's I
        numerator = np.sum(W_sub * np.outer(x_centered, x_centered))
        I = (n_valid / W_sum) * (numerator / sum_sq_dev)

        # Simple permutation test for p-value
        n_perms = 100
        permuted_Is = []
        for _ in range(n_perms):
            perm_scores = np.random.permutation(scores_clean)
            perm_centered = perm_scores - np.mean(perm_scores)
            perm_sum_sq = np.sum(perm_centered**2)
            if perm_sum_sq > 0:
                perm_numerator = np.sum(W_sub * np.outer(perm_centered, perm_centered))
                perm_I = (n_valid / W_sum) * (perm_numerator / perm_sum_sq)
                permuted_Is.append(perm_I)

        if len(permuted_Is) > 0:
            # P-value: proportion of permuted I values >= observed I
            pval = np.mean(np.array(permuted_Is) >= I)
        else:
            pval = np.nan

        moran_results[sig_col] = {"I": I, "pval": pval}

    # Create DataFrame
    moran_df = pd.DataFrame(moran_results).T
    moran_df.columns = ["Morans_I", "p_value"]

    # Create matrix for heatmap (diagonal only)
    moran_I_matrix = pd.DataFrame(index=valid_sigs, columns=valid_sigs, dtype=float)
    for sig in valid_sigs:
        moran_I_matrix.loc[sig, sig] = moran_results[sig]["I"]

    # Generate heatmap
    if not moran_I_matrix.isna().all().all():
        print("\nGenerating Moran's I heatmap...")
        fig, ax = plt.subplots(
            figsize=(max(12, len(valid_sigs) * 0.5), max(10, len(valid_sigs) * 0.5))
        )
        sns.heatmap(
            moran_I_matrix,
            cmap="viridis",
            square=True,
            cbar_kws={"label": "Moran's I"},
            ax=ax,
            fmt=".3f",
            annot=True,
        )
        ax.set_title("Spatial Autocorrelation: Moran's I", fontsize=14, fontweight="bold")
        plt.tight_layout()
        plt.savefig(
            "figures/process_colocalization/morans_i_heatmap.pdf", bbox_inches="tight", dpi=300
        )
        plt.savefig(
            "figures/process_colocalization/morans_i_heatmap.png", bbox_inches="tight", dpi=300
        )
        plt.close()
        print("   ✅ Moran's I heatmap saved")

    # Save results
    moran_df.to_csv(os.path.join(output_dir, "morans_i_results.csv"))
    print(f"   ✅ Results saved to: {output_dir}/morans_i_results.csv")

    # Print summary
    print("\nTop 5 signatures with highest spatial autocorrelation:")
    top_moran = moran_df.nlargest(5, "Morans_I")
    for sig, row in top_moran.iterrows():
        print(f"   {sig}: I = {row['Morans_I']:.3f}, p = {row['p_value']:.3e}")

    return moran_df


def analyze_neighborhood_enrichment(adata, sig_columns, output_dir):
    """
    Analyze thresholded neighborhood enrichment for processes.

    Parameters
    ----------
    adata : AnnData
        Annotated data object with signature scores
    sig_columns : list
        List of signature columns to analyze
    output_dir : str
        Output directory for results

    Returns
    -------
    pd.DataFrame
        DataFrame with neighborhood enrichment results
    """
    print("\n" + "=" * 70)
    print("ANALYSIS 3: NEIGHBORHOOD ENRICHMENT (THRESHOLDED)")
    print("=" * 70)

    # Extract signature scores
    sig_df = adata.obs[sig_columns].copy()
    valid_sigs = sig_df.columns[sig_df.isna().sum() < sig_df.shape[0] * 0.5].tolist()
    print(f"\nAnalyzing {len(valid_sigs)} signatures")
    print(f"   Threshold: low <= {THRESHOLD_LOW}, high >= {THRESHOLD_HIGH}")

    # Ensure spatial neighbors are computed
    if "spatial_neighbors" not in adata.obsp:
        print("\nComputing spatial neighbors...")
        sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)

    # Compute neighborhood enrichment for each signature
    print(f"\nComputing neighborhood enrichment (n_perms={N_PERMS})...")
    nhood_results = {}
    temp_adata = adata.copy()

    for sig in tqdm(valid_sigs, desc="   Computing enrichment"):
        try:
            # Threshold signature
            cat = threshold_signature(temp_adata, sig, THRESHOLD_LOW, THRESHOLD_HIGH)
            cat_binary = (cat == "high").astype(str)

            # Store as categorical
            temp_adata.obs[f"_temp_{sig}"] = pd.Categorical(cat_binary)

            # Compute neighborhood enrichment
            sq.gr.nhood_enrichment(
                temp_adata,
                cluster_key=f"_temp_{sig}",
                n_perms=N_PERMS,
            )

            if "nhood_enrichment" in temp_adata.uns:
                enrichment = temp_adata.uns["nhood_enrichment"]
                if (
                    isinstance(enrichment, pd.DataFrame)
                    and "True" in enrichment.index
                    and "True" in enrichment.columns
                ):
                    nhood_results[sig] = enrichment.loc["True", "True"]
                else:
                    nhood_results[sig] = np.nan
            else:
                nhood_results[sig] = np.nan
        except Exception as e:
            print(f"   Warning: Neighborhood enrichment failed for {sig}: {e}")
            nhood_results[sig] = np.nan

    # Create DataFrame
    nhood_df = pd.DataFrame.from_dict(nhood_results, orient="index", columns=["nhood_enrichment"])

    # Generate heatmap
    if not nhood_df["nhood_enrichment"].isna().all():
        print("\nGenerating neighborhood enrichment heatmap...")
        # Create a matrix (diagonal only for now)
        nhood_matrix = pd.DataFrame(index=valid_sigs, columns=valid_sigs, dtype=float)
        for sig in valid_sigs:
            if sig in nhood_results:
                nhood_matrix.loc[sig, sig] = nhood_results[sig]

        fig, ax = plt.subplots(
            figsize=(max(12, len(valid_sigs) * 0.5), max(10, len(valid_sigs) * 0.5))
        )
        sns.heatmap(
            nhood_matrix,
            cmap="viridis",
            square=True,
            cbar_kws={"label": "Neighborhood Enrichment"},
            ax=ax,
            fmt=".3f",
            annot=True,
        )
        ax.set_title("Neighborhood Enrichment (High vs Not-High)", fontsize=14, fontweight="bold")
        plt.tight_layout()
        plt.savefig(
            "figures/process_colocalization/nhood_enrichment_heatmap.pdf",
            bbox_inches="tight",
            dpi=300,
        )
        plt.savefig(
            "figures/process_colocalization/nhood_enrichment_heatmap.png",
            bbox_inches="tight",
            dpi=300,
        )
        plt.close()
        print("   ✅ Neighborhood enrichment heatmap saved")

    # Save results
    nhood_df.to_csv(os.path.join(output_dir, "neighborhood_enrichment.csv"))
    print(f"   ✅ Results saved to: {output_dir}/neighborhood_enrichment.csv")

    # Print summary
    print("\nTop 5 signatures with highest neighborhood enrichment:")
    top_nhood = nhood_df.nlargest(5, "nhood_enrichment")
    for sig, row in top_nhood.iterrows():
        print(f"   {sig}: enrichment = {row['nhood_enrichment']:.3f}")

    return nhood_df


def analyze_volcano_plots(
    adata, sig_columns, output_dir, signatures_include=None, signatures_exclude=None
):
    """
    Create volcano plots comparing processes across Normal, Non-Solid, and Solid tumor types.

    Parameters
    ----------
    adata : AnnData
        Annotated data object with signature scores
    sig_columns : list
        List of signature columns to analyze
    output_dir : str
        Output directory for results
    signatures_include : list or None, optional
        List of signatures to include. Can be exact names or patterns.
        None means include all. See filter_signatures() for details.
    signatures_exclude : list or None, optional
        List of signatures to exclude. Can be exact names or patterns.
        None means exclude none. See filter_signatures() for details.

    Returns
    -------
    pd.DataFrame
        DataFrame with statistical results
    """
    print("\n" + "=" * 70)
    print("ANALYSIS 4: VOLCANO PLOTS (Normal vs Non-Solid vs Solid)")
    print("=" * 70)

    # Check for tumor_type column
    if "tumor_type" not in adata.obs.columns:
        # Create tumor type mapping
        pa = {
            "Solid": "Solid Tumor",
            "Non-Solid": "Non-Solid Tumor",
            "Normal": "Normal Alveolar Cells",
            "Solid Blood Vessel": "Solid Tumor",
            "Solid Bronchus": "Solid Tumor",
            "Solid Scar Tissue": "Solid Tumor",
            "Non-Solid Blood Vessel": "Non-Solid Tumor",
            "Non-Solid Bronchus": "Non-Solid Tumor",
            "Normal Blood Vessel": "Normal Alveolar Cells",
            "Normal Bronchus": "Normal Alveolar Cells",
            "Solid TLS": "Solid Tumor",
            "TLS Solid": "Solid Tumor",
            "Non-Solid TLS": "Non-Solid Tumor",
            "TLS Non-Solid": "Non-Solid Tumor",
            "TLS Normal": "Normal Alveolar Cells",
        }

        def keep_first_unique(input_list):
            seen = set()
            unique_list = []
            for item in input_list:
                if item not in seen:
                    unique_list.append(item)
                    seen.add(item)
            return unique_list

        adata.obs["tumor_type"] = pd.Categorical(
            adata.obs["pathologist_annotation"].astype(str).replace(pa),
            categories=keep_first_unique(
                ["Normal Alveolar Cells", "Non-Solid Tumor", "Solid Tumor"]
            ),
        )

    # Filter to only the three main categories
    adata_sub = adata[
        adata.obs["tumor_type"].isin(["Normal Alveolar Cells", "Non-Solid Tumor", "Solid Tumor"])
    ].copy()
    adata_sub.obs["tumor_type"] = adata_sub.obs["tumor_type"].cat.remove_unused_categories()

    # Extract signature scores
    sig_df = adata_sub.obs[sig_columns].copy()
    valid_sigs = sig_df.columns[sig_df.isna().sum() < sig_df.shape[0] * 0.5].tolist()

    # Filter signatures based on include/exclude criteria
    original_count = len(valid_sigs)
    if signatures_include is not None or signatures_exclude is not None:
        valid_sigs = filter_signatures(
            valid_sigs, include=signatures_include, exclude=signatures_exclude
        )
        filtered_count = len(valid_sigs)
        print("\nSignature filtering:")
        if signatures_include is not None:
            print(f"   Include: {signatures_include}")
        if signatures_exclude is not None:
            print(f"   Exclude: {signatures_exclude}")
        print(f"   {original_count} -> {filtered_count} signatures")

        if filtered_count == 0:
            print("   ⚠ No signatures remain after filtering!")
            return pd.DataFrame()

        # Show which signatures matched
        print("   Selected signatures:")
        for sig in valid_sigs[:15]:  # Show first 15
            sig_clean = sig.replace("sig:", "").replace("_z", "")
            print(f"      - {sig_clean}")
        if len(valid_sigs) > 15:
            print(f"      ... and {len(valid_sigs) - 15} more")

    print(f"\nAnalyzing {len(valid_sigs)} signatures across 3 tumor types")

    # Define comparisons (reordered to match user request: 1. Normal vs Non-Solid, 2. Non-Solid vs Solid, 3. Normal vs Solid)
    comparisons = [
        ("Normal Alveolar Cells", "Non-Solid Tumor", "Normal_vs_NonSolid", "Normal vs Non-Solid"),
        ("Non-Solid Tumor", "Solid Tumor", "NonSolid_vs_Solid", "Non-Solid vs Solid"),
        ("Normal Alveolar Cells", "Solid Tumor", "Normal_vs_Solid", "Normal vs Solid"),
    ]

    all_results = []
    all_plot_data = []  # Store plot data for faceting

    # First pass: compute all statistics
    for group1, group2, comp_name, comp_title in comparisons:
        print(f"\nComputing statistics for {comp_name}...")
        results = []

        for sig in tqdm(valid_sigs, desc=f"   {comp_name}"):
            # Get values for each group
            group1_values = adata_sub.obs.loc[adata_sub.obs["tumor_type"] == group1, sig].dropna()
            group2_values = adata_sub.obs.loc[adata_sub.obs["tumor_type"] == group2, sig].dropna()

            if len(group1_values) < 3 or len(group2_values) < 3:
                results.append(
                    {
                        "signature": sig,
                        "log2FC": np.nan,
                        "pval": np.nan,
                        "adj_pval": np.nan,
                        "mean_group1": group1_values.mean() if len(group1_values) > 0 else np.nan,
                        "mean_group2": group2_values.mean() if len(group2_values) > 0 else np.nan,
                    }
                )
                continue

            # Compute fold change (log2 of mean ratio)
            mean1 = group1_values.mean()
            mean2 = group2_values.mean()

            # Avoid log(0) or division by zero
            if mean1 == 0 and mean2 == 0:
                log2fc = 0.0
            elif mean1 == 0:
                log2fc = np.sign(mean2) * 10  # Large positive or negative
            elif mean2 == 0:
                log2fc = np.sign(mean1) * -10
            else:
                # Use log2 of ratio, but handle negative values
                if mean1 > 0 and mean2 > 0:
                    log2fc = np.log2(mean2 / mean1)
                else:
                    # For z-scores, use difference instead
                    log2fc = mean2 - mean1

            # Statistical test (Mann-Whitney U)
            try:
                stat, pval = mannwhitneyu(group1_values, group2_values, alternative="two-sided")
            except Exception:
                pval = np.nan

            results.append(
                {
                    "signature": sig,
                    "log2FC": log2fc,
                    "pval": pval,
                    "adj_pval": np.nan,  # Will fill after FDR correction
                    "mean_group1": mean1,
                    "mean_group2": mean2,
                }
            )

        # Convert to DataFrame
        comp_df = pd.DataFrame(results)

        # Apply FDR correction
        valid_pvals = comp_df["pval"].dropna()
        if len(valid_pvals) > 0:
            _, adj_pvals, _, _ = multipletests(valid_pvals, method="fdr_bh", alpha=0.05)
            comp_df.loc[comp_df["pval"].notna(), "adj_pval"] = adj_pvals

        comp_df["comparison"] = comp_name
        comp_df["comp_title"] = comp_title
        comp_df["group1"] = group1
        comp_df["group2"] = group2
        all_results.append(comp_df)

        # Prepare plot data
        plot_df = comp_df.copy()
        plot_df = plot_df.dropna(subset=["log2FC", "pval"])
        plot_df["neg_log10_pval"] = -np.log10(plot_df["pval"] + 1e-300)
        plot_df["color"] = "gray"
        plot_df.loc[(plot_df["adj_pval"] < 0.05) & (plot_df["log2FC"] > 0.5), "color"] = "red"
        plot_df.loc[(plot_df["adj_pval"] < 0.05) & (plot_df["log2FC"] < -0.5), "color"] = "blue"
        plot_df.loc[(plot_df["adj_pval"] < 0.05) & (plot_df["log2FC"].abs() <= 0.5), "color"] = (
            "orange"
        )
        all_plot_data.append(plot_df)

    # Create faceted volcano plot
    print("\n   Generating faceted volcano plots...")
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Plot each comparison in a subplot
    for ax_idx, (comp_name, comp_title, group1, group2) in enumerate(
        [(c[2], c[3], c[0], c[1]) for c in comparisons]
    ):
        ax = axes[ax_idx]
        plot_df = all_plot_data[ax_idx]  # Use the pre-prepared plot data for this comparison

        # Plot
        for color in ["gray", "orange", "red", "blue"]:
            mask = plot_df["color"] == color
            if mask.sum() > 0:
                ax.scatter(
                    plot_df.loc[mask, "log2FC"],
                    plot_df.loc[mask, "neg_log10_pval"],
                    alpha=0.6,
                    s=50,
                    color=color,
                    label={
                        "gray": "Not significant",
                        "orange": "Significant, |FC| ≤ 0.5",
                        "red": "Significant, FC > 0.5",
                        "blue": "Significant, FC < -0.5",
                    }[color]
                    if ax_idx == 0
                    else None,  # Only show legend on first subplot
                )

        # Annotate significant points (only red and blue, not orange)
        # Filter to only red (up) and blue (down) significant points
        annotate_mask = (plot_df["adj_pval"] < 0.05) & (plot_df["color"].isin(["red", "blue"]))
        if annotate_mask.sum() > 0:
            annotate_df = plot_df[annotate_mask].copy()

            # Clean up signature names for display (remove prefix and suffix)
            def clean_sig_name(sig):
                # Remove 'sig:' prefix and '_z' suffix
                cleaned = sig.replace("sig:", "").replace("_z", "")
                # Replace '/' with '-' for readability
                cleaned = cleaned.replace("/", "-")
                # Truncate if too long
                if len(cleaned) > 40:
                    cleaned = cleaned[:37] + "..."
                return cleaned

            annotate_df["label"] = annotate_df["signature"].apply(clean_sig_name)

            # Improved dodging algorithm: stack all labels vertically
            # Sort by y-value (highest p-value first) for consistent ordering
            annotate_df = annotate_df.sort_values("neg_log10_pval", ascending=False).reset_index(
                drop=True
            )

            # Fixed spacing between labels (in points) - ensures no overlap
            label_spacing = 22  # Points between labels (increased for better separation)

            # Group points by x-proximity to stack them vertically
            # Points within this x-distance will be stacked in the same column
            xlim = ax.get_xlim()
            x_range = xlim[1] - xlim[0]
            x_threshold = x_range * 0.10  # 10% of x-range for grouping

            # Group annotations by x-position
            groups = []
            used_indices = set()

            for idx, row in annotate_df.iterrows():
                if idx in used_indices:
                    continue

                x_data = row["log2FC"]
                y_data = row["neg_log10_pval"]

                # Start a new group with this point
                current_group = [
                    {
                        "idx": idx,
                        "x_data": x_data,
                        "y_data": y_data,
                        "label": row["label"],
                        "color": row["color"],
                    }
                ]
                used_indices.add(idx)

                # Find all other points close in x-direction
                for other_idx, other_row in annotate_df.iterrows():
                    if other_idx in used_indices:
                        continue
                    other_x = other_row["log2FC"]
                    if abs(x_data - other_x) < x_threshold:
                        current_group.append(
                            {
                                "idx": other_idx,
                                "x_data": other_x,
                                "y_data": other_row["neg_log10_pval"],
                                "label": other_row["label"],
                                "color": other_row["color"],
                            }
                        )
                        used_indices.add(other_idx)

                # Sort group by y-value (highest first)
                current_group.sort(key=lambda a: a["y_data"], reverse=True)
                groups.append(current_group)

            # For each group, stack labels vertically with consistent spacing
            all_annotations = []
            for group_idx, group in enumerate(groups):
                # Base y-offset for this group
                base_y_offset = 5

                # Stack all labels in this group vertically
                for i, ann in enumerate(group):
                    # Each label gets progressively more y-offset
                    y_offset = base_y_offset + (i * label_spacing)

                    # X-offset: slight variation per group
                    x_offset = 5 + (group_idx % 2) * 3  # Alternate slightly between groups

                    all_annotations.append(
                        {
                            "x_data": ann["x_data"],
                            "y_data": ann["y_data"],
                            "label": ann["label"],
                            "x_offset": x_offset,
                            "y_offset": y_offset,
                            "group_idx": group_idx,
                            "group_pos": i,
                        }
                    )

            # Final pass: ensure all labels are properly spaced
            # Sort by y-data value (highest first)
            all_annotations.sort(key=lambda a: a["y_data"], reverse=True)

            # Adjust y-offsets to ensure minimum spacing
            for i in range(1, len(all_annotations)):
                prev_ann = all_annotations[i - 1]
                curr_ann = all_annotations[i]

                # Calculate what the previous label's final y position would be
                # (approximate: data y + offset in data coordinates)
                # Convert offset from points to data coordinates (rough approximation)
                ylim = ax.get_ylim()
                y_range = ylim[1] - ylim[0]
                fig_ax = ax.figure
                _, fig_height = fig_ax.get_size_inches() * fig_ax.dpi
                points_to_data = y_range / fig_height

                prev_y_final = prev_ann["y_data"] + (prev_ann["y_offset"] * points_to_data)
                curr_y_final = curr_ann["y_data"] + (curr_ann["y_offset"] * points_to_data)

                # Required minimum spacing in data coordinates
                min_spacing_data = label_spacing * points_to_data

                # If labels are too close, increase current y-offset
                if curr_y_final > prev_y_final - min_spacing_data:
                    # Calculate how much more offset is needed
                    needed_spacing = min_spacing_data - (prev_y_final - curr_y_final)
                    needed_offset_points = needed_spacing / points_to_data
                    curr_ann["y_offset"] = prev_ann["y_offset"] + label_spacing

            # Annotate points with dodged positions
            for ann in all_annotations:
                ax.annotate(
                    ann["label"],
                    xy=(ann["x_data"], ann["y_data"]),
                    xytext=(ann["x_offset"], ann["y_offset"]),
                    textcoords="offset points",
                    fontsize=7,
                    alpha=0.8,
                    bbox=dict(
                        boxstyle="round,pad=0.3",
                        facecolor="white",
                        edgecolor="gray",
                        alpha=0.7,
                        linewidth=0.5,
                    ),
                    arrowprops=dict(
                        arrowstyle="->", connectionstyle="arc3,rad=0", lw=0.5, alpha=0.5
                    ),
                )

        # Add significance threshold lines
        ax.axhline(y=-np.log10(0.05), color="black", linestyle="--", linewidth=1, alpha=0.5)
        ax.axvline(x=0.5, color="black", linestyle="--", linewidth=1, alpha=0.3)
        ax.axvline(x=-0.5, color="black", linestyle="--", linewidth=1, alpha=0.3)

        # Labels
        ax.set_xlabel(f"Log2 Fold Change ({group2} / {group1})", fontsize=11)
        ax.set_ylabel("-Log10(p-value)", fontsize=11)
        ax.set_title(comp_title, fontsize=12, fontweight="bold")
        ax.grid(alpha=0.3, linestyle="--")

    # Add legend only on the first subplot
    axes[0].legend(fontsize=8, frameon=True, loc="upper right")

    plt.tight_layout()
    plt.savefig("figures/process_colocalization/volcano_faceted.pdf", bbox_inches="tight", dpi=300)
    plt.savefig("figures/process_colocalization/volcano_faceted.png", bbox_inches="tight", dpi=300)
    plt.close()
    print("   ✅ Faceted volcano plot saved: volcano_faceted.pdf")

    # Combine all results
    all_results_df = pd.concat(all_results, ignore_index=True)

    # Save results
    all_results_df.to_csv(os.path.join(output_dir, "volcano_plot_results.csv"), index=False)
    print(f"\n   ✅ Results saved to: {output_dir}/volcano_plot_results.csv")

    # Print summary
    print("\nSummary of significant differences (adj_p < 0.05, |FC| > 0.5):")
    for comp_name in [c[2] for c in comparisons]:
        comp_data = all_results_df[all_results_df["comparison"] == comp_name]
        sig_up = comp_data[(comp_data["adj_pval"] < 0.05) & (comp_data["log2FC"] > 0.5)]
        sig_down = comp_data[(comp_data["adj_pval"] < 0.05) & (comp_data["log2FC"] < -0.5)]
        print(f"   {comp_name}: {len(sig_up)} up, {len(sig_down)} down")

    return all_results_df


# ============================================================================
# Main Workflow
# ============================================================================


def main():
    """Main function to perform process colocalization analysis."""
    print("=" * 70)
    print("SPOT PROCESS COLOCALIZATION ANALYSIS")
    print("=" * 70)

    # Print configuration
    print("\nAnalysis Configuration:")
    for analysis, enabled in ANALYSIS_CONFIG.items():
        status = "✓ ENABLED" if enabled else "✗ DISABLED"
        print(f"   {analysis}: {status}")

    # Load data
    print("\n1. Loading data...")
    adata_path = Path("results/adata.normalized.scored.p35.h5ad")
    if not adata_path.exists():
        raise FileNotFoundError(f"AnnData file not found: {adata_path}")

    adata = sc.read_h5ad(adata_path)
    print(f"   Loaded: {adata.shape} (spots x genes)")

    # Check for spatial coordinates
    if "spatial" not in adata.obsm:
        raise ValueError("Spatial coordinates not found. Expected adata.obsm['spatial']")

    print(f"   Spatial coordinates: {adata.obsm['spatial'].shape}")

    # Get signature columns
    sig_columns = get_signature_columns(adata)
    print(f"\n2. Found {len(sig_columns)} signature columns")

    if len(sig_columns) == 0:
        print("   ⚠ No signature columns found. Please run score_gene_signatures.py first.")
        return

    # Create output directory
    output_dir = "results/process_colocalization"
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs("figures/process_colocalization", exist_ok=True)

    # Run enabled analyses
    results = {}

    if ANALYSIS_CONFIG["pearson_correlation"]:
        results["correlation"] = analyze_pearson_correlation(adata, sig_columns, output_dir)

    if ANALYSIS_CONFIG["morans_i"]:
        results["morans_i"] = analyze_morans_i(adata, sig_columns, output_dir)

    if ANALYSIS_CONFIG["neighborhood_enrichment"]:
        results["nhood_enrichment"] = analyze_neighborhood_enrichment(
            adata, sig_columns, output_dir
        )

    if ANALYSIS_CONFIG["volcano_plots"]:
        results["volcano"] = analyze_volcano_plots(
            adata,
            sig_columns,
            output_dir,
            signatures_include=VOLCANO_SIGNATURES_INCLUDE,
            signatures_exclude=VOLCANO_SIGNATURES_EXCLUDE,
        )

    # Final summary
    print("\n" + "=" * 70)
    print("✅ PROCESS COLOCALIZATION ANALYSIS COMPLETE")
    print("=" * 70)

    if ANALYSIS_CONFIG["pearson_correlation"]:
        print(f"   ✓ Correlation matrix: {output_dir}/correlation_matrix.csv")
        print("   ✓ Correlation heatmap: figures/process_colocalization/correlation_heatmap.pdf")
        print(
            "   ✓ Correlation clustermap: figures/process_colocalization/correlation_clustermap.pdf"
        )
        print(
            "   ✓ Correlation clustermap: figures/process_colocalization/correlation_clustermap.pdf"
        )

    if ANALYSIS_CONFIG["morans_i"]:
        print(f"   ✓ Moran's I results: {output_dir}/morans_i_results.csv")
        print("   ✓ Moran's I heatmap: figures/process_colocalization/morans_i_heatmap.pdf")

    if ANALYSIS_CONFIG["neighborhood_enrichment"]:
        print(f"   ✓ Neighborhood enrichment: {output_dir}/neighborhood_enrichment.csv")
        print(
            "   ✓ Enrichment heatmap: figures/process_colocalization/nhood_enrichment_heatmap.pdf"
        )

    if ANALYSIS_CONFIG["volcano_plots"]:
        print(f"   ✓ Volcano plot results: {output_dir}/volcano_plot_results.csv")
        print("   ✓ Faceted volcano plot: figures/process_colocalization/volcano_faceted.pdf")

    print("=" * 70)

    return results


if __name__ == "__main__":
    main()
