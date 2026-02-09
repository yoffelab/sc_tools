"""
Macrophage Localization: Measure colocalization with proliferative program.

Tasks:
1. Measure extent of macrophage colocalization with proliferative program
2. Stratify normal, non-solid, and solid tumors
3. Plot tumor proliferation score (x-axis) vs Liron/macrophage scores (y-axis)
4. Scatterplot with different colors per condition (tumor type)
5. Annotate Pearson correlation in separate lines per condition outside of plot
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from scipy import stats
from pathlib import Path
import os

sc.settings.set_figure_params(dpi=300, dpi_save=400)

# ============================================================================
# Configuration
# ============================================================================

# Proliferative program
PROLIFERATIVE_SIG = 'sig:Tumor_Cells/Proliferative_Tumor_z'

# Lung Adenocarcinoma signature (for filtering)
LUNG_ADENOCARCINOMA_SIG = 'sig:Tumor_Cells/Lung Adenocarcinoma_z'

# Macrophage signatures from Liron
MACROPHAGE_SIGNATURES = [
    'sig:Liron/M2-macrophages_z',
    'sig:Liron/MoMacs Merad_z',
    'sig:Liron/Alveolar macrophages Merad_z',
    'sig:Liron/Mac.2 MoMac M2-like_z',
    'sig:Liron/Mac.6 MoMAc M2-like_z',
]

# Tumor type mapping
TUMOR_TYPE_MAPPING = {
    'Solid': 'Solid Tumor',
    'Non-Solid': 'Non-Solid Tumor',
    'Normal': 'Normal Alveolar Cells',
    'Solid Blood Vessel': 'Solid Tumor',
    'Solid Bronchus': 'Solid Tumor',
    'Solid Scar Tissue': 'Solid Tumor',
    'Non-Solid Blood Vessel': 'Non-Solid Tumor',
    'Non-Solid Bronchus': 'Non-Solid Tumor',
    'Normal Blood Vessel': 'Normal Alveolar Cells',
    'Normal Bronchus': 'Normal Alveolar Cells',
    'Solid TLS': 'Solid Tumor',
    'TLS Solid': 'Solid Tumor',
    'Non-Solid TLS': 'Non-Solid Tumor',
    'TLS Non-Solid': 'Non-Solid Tumor',
    'TLS Normal': 'Normal Alveolar Cells',
}

# Colors for tumor types
TUMOR_COLORS = {
    'Normal Alveolar Cells': '#66c2a5',
    'Non-Solid Tumor': '#fc8d62',
    'Solid Tumor': '#8da0cb',
}

# ============================================================================
# Helper Functions
# ============================================================================

def keep_first_unique(input_list):
    """Keep first occurrence of each unique item."""
    seen = set()
    unique_list = []
    for item in input_list:
        if item not in seen:
            unique_list.append(item)
            seen.add(item)
    return unique_list


def calculate_correlation(x, y):
    """Calculate Pearson correlation, handling NaN values."""
    mask = ~(np.isnan(x) | np.isnan(y))
    if np.sum(mask) < 3:  # Need at least 3 points
        return np.nan, np.nan
    r, p = pearsonr(x[mask], y[mask])
    return r, p


def plot_macrophage_colocalization(adata, macrophage_sig, proliferative_sig='sig:Tumor_Cells/Proliferative_Tumor_z',
                                   group_col='tumor_type', output_dir='figures/manuscript/macrophage_localization'):
    """
    Plot scatterplot of proliferative score vs macrophage score, stratified by tumor type.
    Includes filtering and regression lines.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data object with signature scores
    macrophage_sig : str
        Macrophage signature column name
    proliferative_sig : str
        Proliferative signature column name
    group_col : str
        Column name for tumor type grouping
    output_dir : str
        Directory to save plots
    """
    # Check if signatures exist
    if proliferative_sig not in adata.obs.columns:
        print(f"Warning: {proliferative_sig} not found in data. Skipping...")
        return
    
    if macrophage_sig not in adata.obs.columns:
        print(f"Warning: {macrophage_sig} not found in data. Skipping...")
        return
    
    # Filter to only the three main categories
    valid_types = ['Normal Alveolar Cells', 'Non-Solid Tumor', 'Solid Tumor']
    adata_sub = adata[adata.obs[group_col].isin(valid_types)].copy()
    adata_sub.obs[group_col] = adata_sub.obs[group_col].cat.remove_unused_categories()
    
    if adata_sub.n_obs == 0:
        print(f"Warning: No data after filtering. Skipping {macrophage_sig}...")
        return
    
    # Extract data
    x_data = adata_sub.obs[proliferative_sig].values
    y_data = adata_sub.obs[macrophage_sig].values
    groups = adata_sub.obs[group_col].values
    
    # Filter: keep only spots with proliferative score >= 0 AND macrophage score >= 0
    filter_mask = (x_data >= 0) & (y_data >= 0)
    x_filtered = x_data[filter_mask]
    y_filtered = y_data[filter_mask]
    groups_filtered = groups[filter_mask]
    
    # Points to exclude (proliferative score < 0 OR macrophage score < 0)
    x_excluded = x_data[~filter_mask]
    y_excluded = y_data[~filter_mask]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(7, 6))
    
    # Plot excluded points in light gray
    if len(x_excluded) > 0:
        ax.scatter(x_excluded, y_excluded, alpha=0.3, s=15, color='lightgray', 
                  edgecolors='none', label='_nolegend_', zorder=0)
    
    # Add x=0 and y=0 lines
    ax.axvline(x=0, color='black', linestyle='--', linewidth=0.8, alpha=0.5, zorder=1)
    ax.axhline(y=0, color='black', linestyle='--', linewidth=0.8, alpha=0.5, zorder=1)
    
    # Plot scatter and regression for each tumor type
    correlations = {}
    for tumor_type in valid_types:
        if tumor_type not in adata_sub.obs[group_col].cat.categories:
            continue
        
        mask = groups_filtered == tumor_type
        x_vals = x_filtered[mask]
        y_vals = y_filtered[mask]
        
        # Remove NaN values
        valid_mask = ~(np.isnan(x_vals) | np.isnan(y_vals))
        x_vals = x_vals[valid_mask]
        y_vals = y_vals[valid_mask]
        
        if len(x_vals) < 3:
            print(f"Warning: Too few points for {tumor_type} ({len(x_vals)}). Skipping...")
            continue
        
        # Calculate correlation
        r, p = calculate_correlation(x_vals, y_vals)
        correlations[tumor_type] = {'r': r, 'p': p, 'n': len(x_vals)}
        
        # Plot scatter
        color = TUMOR_COLORS.get(tumor_type, 'gray')
        label = tumor_type.replace(' Alveolar Cells', '').replace(' Tumor', '')
        ax.scatter(x_vals, y_vals, alpha=0.5, s=20, color=color, label=label, 
                  edgecolors='none', zorder=3)
        
        # Fit regression line (no confidence intervals)
        if len(x_vals) >= 3:
            # Sort for plotting
            sort_idx = np.argsort(x_vals)
            x_sorted = x_vals[sort_idx]
            y_sorted = y_vals[sort_idx]
            
            # Fit linear regression
            try:
                slope, intercept, r_value, p_value, std_err = stats.linregress(x_sorted, y_sorted)
                
                # Check for valid regression results
                if not (np.isfinite(slope) and np.isfinite(intercept) and np.isfinite(std_err)):
                    print(f"Warning: Invalid regression for {tumor_type}, skipping regression line")
                else:
                    # Generate x values for regression line
                    x_min, x_max = x_sorted.min(), x_sorted.max()
                    if x_max > x_min:  # Ensure there's a range
                        x_line = np.linspace(x_min, x_max, 100)
                        y_line = intercept + slope * x_line
                        
                        # Plot regression line only (no confidence intervals)
                        ax.plot(x_line, y_line, color=color, linewidth=2, alpha=0.8, 
                               linestyle='-', zorder=2, label='_nolegend_')
            except Exception as e:
                print(f"Warning: Error fitting regression for {tumor_type}: {e}")
    
    # Format axes
    proliferative_label = proliferative_sig.replace('sig:', '').replace('_z', '').replace('/', ' ')
    macrophage_label = macrophage_sig.replace('sig:Liron/', '').replace('_z', '').replace('MoMac', 'MoMac ').replace('MoMAc', 'MoMac ')
    
    ax.set_xlabel(f'Proliferative Score (z-scored)', fontsize=12)
    ax.set_ylabel(f'{macrophage_label} Score (z-scored)', fontsize=12)
    ax.set_title(f'Macrophage Colocalization with Proliferation', fontsize=13, fontweight='bold')
    ax.grid(alpha=0.3, linestyle='--', zorder=0)
    ax.legend(title='Tumor Type', fontsize=10, title_fontsize=11, frameon=True, fancybox=True, shadow=True)
    sns.despine(ax=ax)
    
    # Add correlation annotations outside plot
    correlation_text = []
    for tumor_type in valid_types:
        if tumor_type in correlations:
            r = correlations[tumor_type]['r']
            p = correlations[tumor_type]['p']
            n = correlations[tumor_type]['n']
            
            if not np.isnan(r):
                # Format correlation
                if p < 0.001:
                    p_str = 'p<0.001'
                elif p < 0.01:
                    p_str = f'p={p:.3f}'
                elif p < 0.05:
                    p_str = f'p={p:.3f}'
                else:
                    p_str = f'p={p:.3f} (n.s.)'
                
                label = tumor_type.replace(' Alveolar Cells', '').replace(' Tumor', '')
                correlation_text.append(f"{label}: r={r:.3f}, {p_str} (n={n})")
    
    # Add text box with correlations
    if correlation_text:
        correlation_str = "\n".join(correlation_text)
        fig.text(0.5, 0.02, correlation_str,
                ha='center', va='bottom', fontsize=9,
                bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow',
                         alpha=0.8, edgecolor='black', linewidth=1),
                family='monospace')
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0.12, 1, 1])
    
    # Save figure
    os.makedirs(output_dir, exist_ok=True)
    sig_name = macrophage_sig.replace('sig:Liron/', '').replace('_z', '').replace('/', '_')
    output_pdf = os.path.join(output_dir, f'{sig_name}_vs_proliferation.pdf')
    output_png = os.path.join(output_dir, f'{sig_name}_vs_proliferation.png')
    
    plt.savefig(output_pdf, bbox_inches='tight', dpi=300)
    plt.savefig(output_png, bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"  Saved: {sig_name}")
    return correlations


# ============================================================================
# Main Workflow
# ============================================================================

def main():
    """Main function to perform macrophage localization analysis."""
    print("="*70)
    print("Macrophage Localization Analysis")
    print("="*70)
    
    # Load data
    print("\n1. Loading data...")
    adata_path = Path('results/adata.img.genescores.h5ad')
    if not adata_path.exists():
        raise FileNotFoundError(f"AnnData file not found: {adata_path}")
    
    adata = sc.read_h5ad(adata_path)
    print(f"   Loaded: {adata.shape} (spots x genes)")
    
    # Create tumor type column
    print("\n2. Creating tumor type annotations...")
    adata.obs['tumor_type'] = pd.Categorical(
        adata.obs['pathologist_annotation'].astype(str).replace(TUMOR_TYPE_MAPPING),
        categories=keep_first_unique(['Normal Alveolar Cells', 'Non-Solid Tumor', 'Solid Tumor'])
    )
    
    # Filter to only the three main categories
    adata = adata[adata.obs['tumor_type'].isin(['Normal Alveolar Cells', 'Non-Solid Tumor', 'Solid Tumor'])].copy()
    adata.obs['tumor_type'] = adata.obs['tumor_type'].cat.remove_unused_categories()
    
    print(f"   Tumor type distribution:")
    print(adata.obs['tumor_type'].value_counts())
    
    # Check available signatures
    print("\n3. Checking available signatures...")
    available_mac = [sig for sig in MACROPHAGE_SIGNATURES if sig in adata.obs.columns]
    missing_mac = [sig for sig in MACROPHAGE_SIGNATURES if sig not in adata.obs.columns]
    
    if available_mac:
        print(f"   ✓ Found {len(available_mac)} macrophage signatures:")
        for sig in available_mac:
            print(f"     - {sig}")
    
    if missing_mac:
        print(f"   ⚠ Missing {len(missing_mac)} macrophage signatures:")
        for sig in missing_mac:
            print(f"     - {sig}")
    
    if PROLIFERATIVE_SIG not in adata.obs.columns:
        print(f"   ✗ Proliferative signature not found: {PROLIFERATIVE_SIG}")
        print("   Please run score_gene_signatures.py first to generate signature scores.")
        return
    
    print(f"   ✓ Proliferative signature found: {PROLIFERATIVE_SIG}")
    
    # Generate plots for each macrophage signature
    print("\n4. Generating colocalization plots...")
    all_correlations = {}
    
    for macrophage_sig in available_mac:
        print(f"\n   Processing: {macrophage_sig}")
        correlations = plot_macrophage_colocalization(
            adata, macrophage_sig, PROLIFERATIVE_SIG, 
            group_col='tumor_type'
        )
        if correlations:
            all_correlations[macrophage_sig] = correlations
    
    # Save correlation summary
    if all_correlations:
        print("\n5. Saving correlation summary...")
        summary_data = []
        for mac_sig, corrs in all_correlations.items():
            for tumor_type, stats in corrs.items():
                summary_data.append({
                    'Macrophage_Signature': mac_sig.replace('sig:Liron/', '').replace('_z', ''),
                    'Tumor_Type': tumor_type,
                    'Pearson_r': stats['r'],
                    'p_value': stats['p'],
                    'n_spots': stats['n']
                })
        
        summary_df = pd.DataFrame(summary_data)
        summary_path = Path('figures/stats/macrophage_proliferation_correlations.csv')
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_csv(summary_path, index=False)
        print(f"   Saved: {summary_path}")
        
        # Create correlation matrix for visualization
        print("\n6. Generating correlation heatmaps...")
        
        # Pivot data to create correlation matrix (macrophage signatures x tumor types)
        corr_matrix = summary_df.pivot(
            index='Macrophage_Signature',
            columns='Tumor_Type',
            values='Pearson_r'
        )
        
        # Order tumor types as Normal, Non-Solid, Solid
        tumor_order = ['Normal Alveolar Cells', 'Non-Solid Tumor', 'Solid Tumor']
        corr_matrix = corr_matrix[[col for col in tumor_order if col in corr_matrix.columns]]
        
        # Create p-value matrix for annotations
        pval_matrix = summary_df.pivot(
            index='Macrophage_Signature',
            columns='Tumor_Type',
            values='p_value'
        )
        pval_matrix = pval_matrix[[col for col in tumor_order if col in pval_matrix.columns]]
        
        # Create annotation matrix with significance stars
        annot_matrix = corr_matrix.copy()
        for i in range(annot_matrix.shape[0]):
            for j in range(annot_matrix.shape[1]):
                r_val = corr_matrix.iloc[i, j]
                p_val = pval_matrix.iloc[i, j]
                if pd.notna(r_val) and pd.notna(p_val):
                    if p_val < 0.001:
                        stars = '***'
                    elif p_val < 0.01:
                        stars = '**'
                    elif p_val < 0.05:
                        stars = '*'
                    else:
                        stars = ''
                    annot_matrix.iloc[i, j] = f"{r_val:.3f}{stars}"
                else:
                    annot_matrix.iloc[i, j] = ""
        
        # Generate regular heatmap
        fig, ax = plt.subplots(figsize=(8, 6))
        sns.heatmap(
            corr_matrix,
            annot=annot_matrix,
            fmt='',
            cmap='RdBu_r',
            center=0,
            vmin=-1,
            vmax=1,
            cbar_kws={'label': 'Pearson r'},
            linewidths=0.5,
            linecolor='gray',
            ax=ax
        )
        ax.set_title('Macrophage-Proliferation Correlations', fontsize=14, fontweight='bold', pad=15)
        ax.set_xlabel('Tumor Type', fontsize=12)
        ax.set_ylabel('Macrophage Signature', fontsize=12)
        
        # Rotate labels for better readability
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
        
        plt.tight_layout()
        heatmap_pdf = Path('figures/manuscript/macrophage_localization/correlation_heatmap.pdf')
        heatmap_png = Path('figures/manuscript/macrophage_localization/correlation_heatmap.png')
        plt.savefig(heatmap_pdf, bbox_inches='tight', dpi=300)
        plt.savefig(heatmap_png, bbox_inches='tight', dpi=300)
        plt.close()
        print(f"   Saved heatmap: {heatmap_pdf.name}")
        
        # Generate clustermap (hierarchically clustered heatmap)
        g = sns.clustermap(
            corr_matrix,
            cmap='RdBu_r',
            center=0,
            vmin=-1,
            vmax=1,
            figsize=(9, 7),
            cbar_kws={'label': 'Pearson r'},
            linewidths=0.5,
            linecolor='gray',
            dendrogram_ratio=0.15,
            cbar_pos=(0.02, 0.8, 0.03, 0.15),
            row_cluster=True,
            col_cluster=False,  # Don't cluster tumor types (keep order: Normal, Non-Solid, Solid)
            yticklabels=True,
            xticklabels=True,
            annot=annot_matrix,
            fmt=''
        )
        
        # Set title
        g.fig.suptitle('Macrophage-Proliferation Correlations (Clustered)', 
                      fontsize=14, fontweight='bold', y=0.98)
        
        # Adjust labels
        g.ax_heatmap.set_xlabel('Tumor Type', fontsize=12)
        g.ax_heatmap.set_ylabel('Macrophage Signature', fontsize=12)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        
        clustermap_pdf = Path('figures/manuscript/macrophage_localization/correlation_clustermap.pdf')
        clustermap_png = Path('figures/manuscript/macrophage_localization/correlation_clustermap.png')
        g.savefig(clustermap_pdf, bbox_inches='tight', dpi=300)
        g.savefig(clustermap_png, bbox_inches='tight', dpi=300)
        plt.close()
        print(f"   Saved clustermap: {clustermap_pdf.name}")
    
    print("\n" + "="*70)
    print("✅ Macrophage localization analysis complete!")
    print(f"   - Scatterplots saved to: figures/manuscript/macrophage_localization/")
    if all_correlations:
        print(f"   - Correlations CSV: figures/stats/macrophage_proliferation_correlations.csv")
        print(f"   - Correlation heatmap: figures/manuscript/macrophage_localization/correlation_heatmap.pdf")
        print(f"   - Correlation clustermap: figures/manuscript/macrophage_localization/correlation_clustermap.pdf")
    print("="*70)


if __name__ == "__main__":
    main()

