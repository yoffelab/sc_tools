"""
TLS-Specific Transcriptomics Analysis.

Tasks:
1. TLS Niche Extraction: Subset tls_clustered.h5ad to focus on lymphoid-rich neighborhoods
2. B-cell/T-cell State Comparison: Analyze ratios of Naive vs. Memory B-cells and 
   Exhausted vs. Effector T-cells within TLS across tumor stages
3. Ligand-Receptor Crosstalk: Infer signaling between tumor-adjacent TLS and tumor cells
   in GGO vs. Solid samples using squidpy.gr.ligrec
"""

import scanpy as sc
import json
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import seaborn as sns
import squidpy as sq
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

sc.settings.set_figure_params(dpi=300, dpi_save=400)

# Load data
print("Loading data...")
with open('metadata/gene_signatures.json', 'r') as file:
    spatial_signatures = json.load(file)

adata = sc.read('results/adata.normalized.scored.p35.h5ad')
adata_tls = sc.read('results/tls_clustered.h5ad')

adata.obs['solidity_type'] = pd.Categorical(
    adata.obs['solidity_type'], 
    categories=['Normal', 'Non-Solid', 'Solid']
)
adata_tls.obs['solidity_type_majority'] = pd.Categorical(
    adata_tls.obs['solidity_type_majority'], 
    categories=['Normal', 'Non-Solid', 'Solid']
)

# ============================================================================
# 1. Basic TLS Visualization
# ============================================================================
print("\n" + "="*70)
print("1. BASIC TLS VISUALIZATION")
print("="*70)

# Plot UMAP of TLSs
print("\nComputing UMAP for TLS clusters...")
if 'X_pca' not in adata_tls.obsm:
    sc.pp.pca(adata_tls)
if 'neighbors' not in adata_tls.obsp:
    sc.pp.neighbors(adata_tls)
sc.tl.umap(adata_tls)

sc.pl.umap(
    adata_tls, 
    color=['library_id', 'solidity_type_majority', 'tls_size_mm2'], 
    show=False, 
    save='tls_umap'
)

# Plot spatial of TLSs
os.makedirs('figures/tls_spatial', exist_ok=True)
print("\nPlotting spatial distribution of TLSs...")
for lib_id in adata.obs['library_id'].unique():
    adata_lib = adata[adata.obs['library_id'] == lib_id].copy()
    
    sc.pl.spatial(
        adata_lib,
        color=['architecture_type', 'solidity_type'], 
        library_id=lib_id,
        show=False, 
        save=None
    )
    plt.savefig(f'figures/tls_spatial/{lib_id}.pdf', bbox_inches='tight')
    plt.savefig(f'figures/tls_spatial/{lib_id}.png', bbox_inches='tight')
    plt.close()

# ============================================================================
# 2. TLS Niche Extraction: Score B-cell and T-cell signatures
# ============================================================================
print("\n" + "="*70)
print("2. TLS NICHE EXTRACTION: B-CELL AND T-CELL STATE SCORING")
print("="*70)

# Function to score gene signatures
def score_signature(adata, genes, use_raw=True):
    """Score a gene signature by averaging z-scored expression."""
    if use_raw and adata.raw is not None:
        X = adata.raw.X
        var_names = adata.raw.var_names
    else:
        X = adata.X
        var_names = adata.var_names
    
    # Find genes present in data
    present_genes = [g for g in genes if g in var_names]
    if len(present_genes) == 0:
        return np.zeros(adata.n_obs)
    
    # Get indices
    gene_indices = [list(var_names).index(g) for g in present_genes]
    
    # Extract expression
    if hasattr(X, 'toarray'):
        X_dense = X.toarray()
    else:
        X_dense = X
    
    # Z-score per gene, then average
    scores = np.zeros(adata.n_obs)
    for idx in gene_indices:
        gene_expr = X_dense[:, idx]
        gene_expr_z = (gene_expr - np.mean(gene_expr)) / (np.std(gene_expr) + 1e-8)
        scores += gene_expr_z
    
    return scores / len(present_genes)

# Score B-cell states
print("\nScoring B-cell states...")
b_naive_genes = spatial_signatures['Immune_Lymphoid']['B_Naive']
b_memory_genes = spatial_signatures['Immune_Lymphoid']['B_Memory']

adata_tls.obs['B_Naive_score'] = score_signature(adata_tls, b_naive_genes, use_raw=True)
adata_tls.obs['B_Memory_score'] = score_signature(adata_tls, b_memory_genes, use_raw=True)

# Calculate ratio (avoid division by zero)
b_total = adata_tls.obs['B_Naive_score'] + adata_tls.obs['B_Memory_score']
adata_tls.obs['B_Naive_ratio'] = np.where(
    b_total > 0,
    adata_tls.obs['B_Naive_score'] / b_total,
    np.nan
)
adata_tls.obs['B_Memory_ratio'] = np.where(
    b_total > 0,
    adata_tls.obs['B_Memory_score'] / b_total,
    np.nan
)

# Score T-cell states
print("Scoring T-cell states...")
t_exhausted_genes = spatial_signatures['Immune_Lymphoid']['Exhausted_CD8_T']
t_effector_genes = spatial_signatures['Immune_Lymphoid']['Effector_CD8_T']

adata_tls.obs['T_Exhausted_score'] = score_signature(adata_tls, t_exhausted_genes, use_raw=True)
adata_tls.obs['T_Effector_score'] = score_signature(adata_tls, t_effector_genes, use_raw=True)

# Calculate ratio
t_total = adata_tls.obs['T_Exhausted_score'] + adata_tls.obs['T_Effector_score']
adata_tls.obs['T_Exhausted_ratio'] = np.where(
    t_total > 0,
    adata_tls.obs['T_Exhausted_score'] / t_total,
    np.nan
)
adata_tls.obs['T_Effector_ratio'] = np.where(
    t_total > 0,
    adata_tls.obs['T_Effector_score'] / t_total,
    np.nan
)

# ============================================================================
# 3. B-cell/T-cell State Comparison Across Tumor Stages
# ============================================================================
print("\n" + "="*70)
print("3. B-CELL/T-CELL STATE COMPARISON ACROSS TUMOR STAGES")
print("="*70)

os.makedirs('figures/tls_spatial/cell_states', exist_ok=True)
os.makedirs('results/tls', exist_ok=True)

# Statistical testing function
def test_one_vs_all(adata, score_col, group_col='solidity_type_majority'):
    """Perform 1-vs-all comparisons with FDR correction."""
    results = []
    groups = adata.obs[group_col].cat.categories
    
    for group in groups:
        group_mask = adata.obs[group_col] == group
        rest_mask = ~group_mask
        
        group_values = adata.obs.loc[group_mask, score_col].dropna()
        rest_values = adata.obs.loc[rest_mask, score_col].dropna()
        
        if len(group_values) < 3 or len(rest_values) < 3:
            results.append({
                'group': group,
                'pval': np.nan,
                'adj_pval': np.nan,
                'n_group': len(group_values),
                'n_rest': len(rest_values),
            })
            continue
        
        stat, pval = mannwhitneyu(group_values, rest_values, alternative='two-sided')
        results.append({
            'group': group,
            'pval': pval,
            'adj_pval': np.nan,
            'n_group': len(group_values),
            'n_rest': len(rest_values),
        })
    
    df = pd.DataFrame(results)
    
    # Apply FDR correction
    valid_pvals = df['pval'].dropna()
    if len(valid_pvals) > 0:
        _, adj_pvals, _, _ = multipletests(valid_pvals, method='fdr_bh', alpha=0.05)
        df.loc[df['pval'].notna(), 'adj_pval'] = adj_pvals
    
    return df

# Test B-cell ratios
print("\nTesting B-cell ratios across tumor stages...")
b_naive_stats = test_one_vs_all(adata_tls, 'B_Naive_ratio', 'solidity_type_majority')
b_memory_stats = test_one_vs_all(adata_tls, 'B_Memory_ratio', 'solidity_type_majority')

# Test T-cell ratios
print("Testing T-cell ratios across tumor stages...")
t_exhausted_stats = test_one_vs_all(adata_tls, 'T_Exhausted_ratio', 'solidity_type_majority')
t_effector_stats = test_one_vs_all(adata_tls, 'T_Effector_ratio', 'solidity_type_majority')

# Combine statistics
all_stats = pd.DataFrame({
    'B_Naive': b_naive_stats.set_index('group')['adj_pval'],
    'B_Memory': b_memory_stats.set_index('group')['adj_pval'],
    'T_Exhausted': t_exhausted_stats.set_index('group')['adj_pval'],
    'T_Effector': t_effector_stats.set_index('group')['adj_pval'],
})
all_stats.to_csv('results/tls/cell_state_ratios_stats.csv')
print(f"\nSaved statistics to results/tls/cell_state_ratios_stats.csv")

# Plot boxplots
def plot_ratio_boxplot(adata, ratio_col, title, stats_df, output_path, figsize=(6, 4)):
    """Plot boxplot with significance annotations."""
    fig, ax = plt.subplots(figsize=figsize)
    
    groups = adata.obs['solidity_type_majority'].cat.categories
    plot_data = [adata.obs.loc[adata.obs['solidity_type_majority'] == g, ratio_col].dropna().values 
                 for g in groups]
    
    bp = ax.boxplot(plot_data, labels=groups, patch_artist=True, showfliers=False, widths=0.6)
    colors = ['#66c2a5', '#fc8d62', '#8da0cb']
    for patch, color in zip(bp['boxes'], colors[:len(groups)]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    # Add significance annotations
    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min
    bar_y = y_max - y_range * 0.05
    
    if ratio_col in stats_df.columns:
        pvals = stats_df[ratio_col].dropna()
        if len(pvals) > 0:
            min_pval = pvals.min()
            ax.text(0.98, 0.98, f'min adj. p = {min_pval:.2e}', 
                   transform=ax.transAxes, ha='right', va='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
                   fontsize=9)
    
    ax.set_ylabel('Ratio', fontsize=11)
    ax.set_xlabel('Tumor Type', fontsize=11)
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    sns.despine(ax=ax)
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()

# Plot B-cell ratios
plot_ratio_boxplot(
    adata_tls, 'B_Naive_ratio', 'B-cell: Naive Ratio',
    all_stats, 'figures/tls_spatial/cell_states/B_Naive_ratio_boxplot.pdf'
)
plot_ratio_boxplot(
    adata_tls, 'B_Memory_ratio', 'B-cell: Memory Ratio',
    all_stats, 'figures/tls_spatial/cell_states/B_Memory_ratio_boxplot.pdf'
)

# Plot T-cell ratios
plot_ratio_boxplot(
    adata_tls, 'T_Exhausted_ratio', 'T-cell: Exhausted Ratio',
    all_stats, 'figures/tls_spatial/cell_states/T_Exhausted_ratio_boxplot.pdf'
)
plot_ratio_boxplot(
    adata_tls, 'T_Effector_ratio', 'T-cell: Effector Ratio',
    all_stats, 'figures/tls_spatial/cell_states/T_Effector_ratio_boxplot.pdf'
)

# ============================================================================
# 4. Ligand-Receptor Crosstalk: Tumor-adjacent TLS vs Tumor cells
# ============================================================================
print("\n" + "="*70)
print("4. LIGAND-RECEPTOR CROSSTALK: TLS-TUMOR INTERACTIONS")
print("="*70)

os.makedirs('results/tls/ligrec', exist_ok=True)
os.makedirs('figures/tls_spatial/ligrec', exist_ok=True)

# Identify tumor-adjacent TLS and tumor cells in spatial data
print("\nIdentifying tumor-adjacent TLS and tumor regions...")

# Create annotation for ligand-receptor analysis
# Mark TLS and tumor regions
adata.obs['region_type'] = 'Other'
tls_mask = adata.obs['architecture_type'].astype(str).str.contains('TLS', case=False, na=False)
tumor_mask = (
    adata.obs['pathologist_annotation'].astype(str).str.contains('Tumor', case=False, na=False) &
    ~tls_mask
)

adata.obs.loc[tls_mask, 'region_type'] = 'TLS'
adata.obs.loc[tumor_mask, 'region_type'] = 'Tumor'
adata.obs['region_type'] = pd.Categorical(adata.obs['region_type'], categories=['TLS', 'Tumor', 'Other'])

# Filter to GGO (Non-Solid) and Solid samples for comparison
for tumor_stage in ['Non-Solid', 'Solid']:
    print(f"\nAnalyzing {tumor_stage} samples...")
    
    stage_mask = adata.obs['solidity_type'] == tumor_stage
    adata_stage = adata[stage_mask].copy()
    
    # Filter to TLS and Tumor regions only
    region_mask = adata_stage.obs['region_type'].isin(['TLS', 'Tumor'])
    if region_mask.sum() < 10:
        print(f"  Skipping {tumor_stage}: insufficient TLS/Tumor regions")
        continue
    
    adata_stage_regions = adata_stage[region_mask].copy()
    
    # Ensure spatial neighbors are computed
    if 'spatial_neighbors' not in adata_stage_regions.obsp:
        print(f"  Computing spatial neighbors for {tumor_stage}...")
        sq.gr.spatial_neighbors(adata_stage_regions, coord_type="generic", delaunay=True)
    
    # Normalize if needed
    if adata_stage_regions.raw is not None:
        adata_stage_regions.X = adata_stage_regions.raw.X
    
    if np.max(adata_stage_regions.X) > 100:
        sc.pp.normalize_total(adata_stage_regions, target_sum=1e4)
        sc.pp.log1p(adata_stage_regions)
    
    # Run ligand-receptor analysis
    try:
        print(f"  Running ligand-receptor analysis for {tumor_stage}...")
        sq.gr.ligrec(  # Note: may be sq.gr.ligprec in some squidpy versions
            adata_stage_regions,
            n_perms=1000,
            cluster_key='region_type',
            copy=True,
        )
        
        if 'ligrec' in adata_stage_regions.uns:
            ligrec_results = adata_stage_regions.uns['ligrec']
            if isinstance(ligrec_results, dict) and 'pvalues' in ligrec_results:
                ligrec_df = ligrec_results['pvalues']
                ligrec_df.to_csv(f'results/tls/ligrec/ligrec_{tumor_stage}.csv')
                print(f"  Saved ligand-receptor results for {tumor_stage}")
        
        # Plot spatial distribution
        for lib_id in adata_stage_regions.obs['library_id'].unique()[:2]:
            adata_lib = adata_stage_regions[adata_stage_regions.obs['library_id'] == lib_id].copy()
            sc.pl.spatial(
                adata_lib,
                color='region_type',
                library_id=lib_id,
                show=False,
                title=f'{lib_id} ({tumor_stage}): TLS-Tumor Regions'
            )
            plt.savefig(
                f'figures/tls_spatial/ligrec/{lib_id}_{tumor_stage}_regions.pdf',
                bbox_inches='tight', dpi=300
            )
            plt.savefig(
                f'figures/tls_spatial/ligrec/{lib_id}_{tumor_stage}_regions.png',
                bbox_inches='tight', dpi=300
            )
            plt.close()
            
    except Exception as e:
        print(f"  ⚠️  Ligand-receptor analysis failed for {tumor_stage}: {e}")

print("\n" + "="*70)
print("✅ TLS-SPECIFIC TRANSCRIPTOMICS ANALYSIS COMPLETE")
print("="*70)
print(f"   - Cell state statistics: results/tls/cell_state_ratios_stats.csv")
print(f"   - Cell state plots: figures/tls_spatial/cell_states/")
print(f"   - Ligand-receptor results: results/tls/ligrec/")
print(f"   - Ligand-receptor plots: figures/tls_spatial/ligrec/")
