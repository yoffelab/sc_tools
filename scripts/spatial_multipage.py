"""
Spatial Multipage PDF: Generate multipage PDF with spatial plots per library.

For each library_id, creates a page with 6 plots (3x2 grid):
1. Plain H&E tissue
2. Pathologist annotation
3. Solidity
4. Proliferation score
5. Macrophage score
6. Truncated similarity (proliferation * macrophage, if both positive, else 0)
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path
import os
from tqdm import tqdm
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import seaborn as sns

sc.settings.set_figure_params(dpi=300, dpi_save=400)

# ============================================================================
# Configuration
# ============================================================================

# Input file
ADATA_PATH = 'results/adata.img.genescores.h5ad'

# Output directory
OUTPUT_DIR = 'figures/manuscript/spatial_multipage'
OUTPUT_PDF = os.path.join(OUTPUT_DIR, 'spatial_multipage.pdf')

# Signature columns
PROLIFERATION_SIG = 'sig:Tumor_Cells/Proliferative_Tumor_z'
MACROPHAGE_SIG = 'sig:Liron/M2-macrophages_z'  # Default, can be changed

# Tumor type mapping (for solidity)
TUMOR_TYPE_MAPPING = {
    'Solid': 'Solid',
    'Non-Solid': 'Non-Solid',
    'Normal': 'Normal',
    'Solid Blood Vessel': 'Solid',
    'Solid Bronchus': 'Solid',
    'Solid Scar Tissue': 'Solid',
    'Non-Solid Blood Vessel': 'Non-Solid',
    'Non-Solid Bronchus': 'Non-Solid',
    'Normal Blood Vessel': 'Normal',
    'Normal Bronchus': 'Normal',
    'Solid TLS': 'Solid',
    'TLS Solid': 'Solid',
    'Non-Solid TLS': 'Non-Solid',
    'TLS Non-Solid': 'Non-Solid',
    'TLS Normal': 'Normal',
}

# Solidity colors
SOLIDITY_COLORS = {
    'Normal': '#66c2a5',
    'Non-Solid': '#fc8d62',
    'Solid': '#8da0cb',
}

# ============================================================================
# Helper Functions
# ============================================================================

def prepare_solidity(adata):
    """Prepare solidity column from pathologist_annotation."""
    if 'solidity_type' in adata.obs.columns:
        return adata.obs['solidity_type']
    
    if 'pathologist_annotation' not in adata.obs.columns:
        raise KeyError("Neither 'solidity_type' nor 'pathologist_annotation' found in adata.obs")
    
    # Map pathologist_annotation to solidity
    solidity = adata.obs['pathologist_annotation'].astype(str).replace(TUMOR_TYPE_MAPPING)
    return pd.Categorical(solidity, categories=['Normal', 'Non-Solid', 'Solid'])


def calculate_truncated_similarity(proliferation, macrophage):
    """
    Calculate truncated similarity: proliferation * macrophage if both positive, else 0.
    
    Parameters
    ----------
    proliferation : array-like
        Proliferation scores
    macrophage : array-like
        Macrophage scores
    
    Returns
    -------
    array-like
        Truncated similarity scores
    """
    proliferation = np.asarray(proliferation)
    macrophage = np.asarray(macrophage)
    
    # Create mask: both positive
    mask = (proliferation > 0) & (macrophage > 0)
    
    # Calculate similarity
    similarity = np.zeros_like(proliferation)
    similarity[mask] = proliferation[mask] * macrophage[mask]
    
    return similarity


def plot_spatial_plain_he(adata_full, adata_sub, library_id, ax):
    """Plot plain H&E tissue image."""
    try:
        # Check if spatial data exists for this library (use full adata, not subset)
        if library_id not in adata_full.uns.get('spatial', {}):
            ax.text(0.5, 0.5, f'No spatial data for library {library_id}', 
                    ha='center', va='center', transform=ax.transAxes)
            ax.set_title('H&E Tissue', fontsize=12, fontweight='bold')
            return
        
        # Get the image
        spatial_data = adata_full.uns['spatial'][library_id]
        if 'images' not in spatial_data or 'hires' not in spatial_data['images']:
            ax.text(0.5, 0.5, f'No H&E image found for library {library_id}', 
                    ha='center', va='center', transform=ax.transAxes)
            ax.set_title('H&E Tissue', fontsize=12, fontweight='bold')
            return
        
        # Display the image
        img = spatial_data['images']['hires']
        ax.imshow(img, aspect='auto')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title('H&E Tissue', fontsize=12, fontweight='bold')
    except Exception as e:
        ax.text(0.5, 0.5, f'Error loading H&E image:\n{str(e)}', 
                ha='center', va='center', transform=ax.transAxes, fontsize=10)
        ax.set_title('H&E Tissue', fontsize=12, fontweight='bold')


def plot_spatial_annotation(adata_sub, library_id, ax):
    """Plot pathologist annotation."""
    if 'pathologist_annotation' not in adata_sub.obs.columns:
        ax.text(0.5, 0.5, 'pathologist_annotation not found', 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Pathologist Annotation', fontsize=12, fontweight='bold')
        return
    
    sc.pl.spatial(
        adata_sub,
        color='pathologist_annotation',
        library_id=library_id,
        frameon=False,
        show=False,
        ax=ax,
        legend_loc='right margin',
    )
    ax.set_title('Pathologist Annotation', fontsize=12, fontweight='bold')


def plot_spatial_solidity(adata_sub, library_id, ax):
    """Plot solidity."""
    if 'solidity' not in adata_sub.obs.columns:
        ax.text(0.5, 0.5, 'solidity not found', 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Solidity', fontsize=12, fontweight='bold')
        return
    
    # Create color map
    unique_solidities = adata_sub.obs['solidity'].cat.categories
    colors = [SOLIDITY_COLORS.get(s, '#808080') for s in unique_solidities]
    
    sc.pl.spatial(
        adata_sub,
        color='solidity',
        library_id=library_id,
        frameon=False,
        show=False,
        ax=ax,
        palette=colors,
        legend_loc='right margin',
    )
    ax.set_title('Solidity', fontsize=12, fontweight='bold')


def plot_spatial_score(adata_sub, library_id, score_col, ax, title, cmap='coolwarm', vmin=None, vmax=None):
    """Plot spatial score."""
    if score_col not in adata_sub.obs.columns:
        ax.text(0.5, 0.5, f'{score_col} not found', 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title, fontsize=12, fontweight='bold')
        return
    
    sc.pl.spatial(
        adata_sub,
        color=score_col,
        library_id=library_id,
        frameon=False,
        show=False,
        ax=ax,
        cmap=cmap,
        colorbar_loc='right',
        vmin=vmin,
        vmax=vmax,
    )
    # Clean up title
    clean_title = title.replace('_', ' ').title()
    ax.set_title(clean_title, fontsize=12, fontweight='bold')


def create_multipage_pdf(adata, output_pdf, proliferation_sig=PROLIFERATION_SIG, 
                         macrophage_sig=MACROPHAGE_SIG):
    """
    Create multipage PDF with spatial plots for each library.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with spatial data (must have solidity and truncated_similarity columns)
    output_pdf : str
        Path to output PDF file
    proliferation_sig : str
        Column name for proliferation score
    macrophage_sig : str
        Column name for macrophage score
    """
    
    # Get unique library IDs
    library_ids = sorted(adata.obs['library_id'].unique())
    print(f"Found {len(library_ids)} libraries: {library_ids}")
    
    # Create output directory
    os.makedirs(os.path.dirname(output_pdf), exist_ok=True)
    
    # Create PDF
    with PdfPages(output_pdf) as pdf:
        for lib_id in tqdm(library_ids, desc="Creating pages"):
            # Subset data for this library
            adata_sub = adata[adata.obs['library_id'] == lib_id].copy()
            
            if len(adata_sub) == 0:
                print(f"Warning: No data for library {lib_id}, skipping...")
                continue
            
            # Create figure with 2x3 grid
            fig, axes = plt.subplots(2, 3, figsize=(18, 12))
            axes = axes.flatten()
            
            # 1. Plain H&E tissue
            plot_spatial_plain_he(adata, adata_sub, lib_id, axes[0])
            
            # 2. Pathologist annotation
            plot_spatial_annotation(adata_sub, lib_id, axes[1])
            
            # 3. Solidity
            plot_spatial_solidity(adata_sub, lib_id, axes[2])
            
            # 4. Proliferation score
            plot_spatial_score(adata_sub, lib_id, proliferation_sig, axes[3], 
                              'Tumor Proliferation', cmap='coolwarm', vmin=-2, vmax=2)
            
            # 5. Macrophage score
            plot_spatial_score(adata_sub, lib_id, macrophage_sig, axes[4], 
                              'M2 Macrophage', cmap='coolwarm', vmin=-2, vmax=2)
            
            # 6. Truncated similarity
            plot_spatial_score(adata_sub, lib_id, 'truncated_similarity', axes[5], 
                              'Truncated Similarity', cmap='YlOrRd', vmin=0, vmax=1)
            
            # Add library ID as suptitle
            fig.suptitle(f'Library: {lib_id}', fontsize=16, fontweight='bold', y=0.995)
            
            # Save page to PDF
            pdf.savefig(fig, bbox_inches='tight', dpi=300)
            plt.close(fig)
    
    print(f"\n✅ Multipage PDF created: {output_pdf}")
    print(f"   Total pages: {len(library_ids)}")


def calculate_proportion_nonzero(adata):
    """
    Calculate proportion of non-zero truncated similarity scores per patient per solidity.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with truncated_similarity and solidity columns
    
    Returns
    -------
    pd.DataFrame
        DataFrame with columns: patient (library_id), solidity, proportion_nonzero
    """
    if 'truncated_similarity' not in adata.obs.columns:
        raise KeyError("truncated_similarity column not found in adata.obs")
    
    if 'solidity' not in adata.obs.columns:
        raise KeyError("solidity column not found in adata.obs")
    
    if 'library_id' not in adata.obs.columns:
        raise KeyError("library_id column not found in adata.obs")
    
    results = []
    
    # Group by patient (library_id) and solidity
    for (patient, solidity), group in adata.obs.groupby(['library_id', 'solidity']):
        if len(group) == 0:
            continue
        
        similarity_scores = group['truncated_similarity'].values
        n_total = len(similarity_scores)
        n_nonzero = np.sum(similarity_scores > 0)
        proportion = n_nonzero / n_total if n_total > 0 else 0.0
        
        results.append({
            'patient': patient,
            'solidity': solidity,
            'proportion_nonzero': proportion,
            'n_spots': n_total,
            'n_nonzero': n_nonzero
        })
    
    return pd.DataFrame(results)


def plot_proportion_boxplot(adata, output_dir):
    """
    Create boxplot of proportion of non-zero truncated similarity per patient per solidity.
    Includes statistical comparisons using Mann-Whitney U test.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with truncated_similarity, solidity, and library_id columns
    output_dir : str
        Output directory for saving plots
    """
    print("\n" + "="*70)
    print("Creating proportion boxplot...")
    print("="*70)
    
    # Calculate proportions
    prop_df = calculate_proportion_nonzero(adata)
    
    if len(prop_df) == 0:
        print("Warning: No data for proportion calculation")
        return
    
    # Filter to only Normal, Non-Solid, Solid
    solidity_order = ['Normal', 'Non-Solid', 'Solid']
    prop_df = prop_df[prop_df['solidity'].isin(solidity_order)].copy()
    prop_df['solidity'] = pd.Categorical(prop_df['solidity'], categories=solidity_order)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Create boxplot with individual points
    sns.boxplot(data=prop_df, x='solidity', y='proportion_nonzero', ax=ax, 
                palette=[SOLIDITY_COLORS.get(s, '#808080') for s in solidity_order],
                width=0.6)
    
    # Add individual points (one per patient)
    sns.stripplot(data=prop_df, x='solidity', y='proportion_nonzero', ax=ax,
                  color='black', size=5, alpha=0.6, jitter=True)
    
    ax.set_xlabel('Solidity', fontsize=12, fontweight='bold')
    ax.set_ylabel('Proportion of Non-Zero Truncated Similarity', fontsize=12, fontweight='bold')
    ax.set_title('Proportion of Non-Zero Truncated Similarity\nby Solidity Condition', 
                 fontsize=14, fontweight='bold', pad=15)
    
    # Perform statistical comparisons
    print("\nPerforming statistical comparisons (Mann-Whitney U test)...")
    comparisons = []
    p_values = []
    
    solidity_groups = {}
    for solidity in solidity_order:
        solidity_groups[solidity] = prop_df[prop_df['solidity'] == solidity]['proportion_nonzero'].values
    
    # All pairwise comparisons
    for i, group1 in enumerate(solidity_order):
        for j, group2 in enumerate(solidity_order):
            if i >= j:
                continue
            
            values1 = solidity_groups[group1]
            values2 = solidity_groups[group2]
            
            if len(values1) < 2 or len(values2) < 2:
                print(f"   ⚠ Skipping {group1} vs {group2}: insufficient data")
                continue
            
            stat, pval = mannwhitneyu(values1, values2, alternative='two-sided')
            comparisons.append(f"{group1} vs {group2}")
            p_values.append(pval)
            print(f"   {group1} vs {group2}: p = {pval:.4f}")
    
    # Apply FDR correction
    if len(p_values) > 0:
        _, p_adjusted, _, _ = multipletests(p_values, method='fdr_bh', alpha=0.05)
        
        # Add significance annotations
        y_max = prop_df['proportion_nonzero'].max()
        y_range = prop_df['proportion_nonzero'].max() - prop_df['proportion_nonzero'].min()
        y_offset = y_max + 0.05 * y_range
        
        # Function to get asterisks
        def get_asterisks(p):
            if p < 0.0001:
                return '****'
            elif p < 0.001:
                return '***'
            elif p < 0.01:
                return '**'
            elif p < 0.05:
                return '*'
            else:
                return 'ns'
        
        # Draw significance bars
        bar_height = 0.02 * y_range
        bar_y = y_offset
        bar_spacing = 0.03 * y_range
        
        for idx, (comp, p_adj) in enumerate(zip(comparisons, p_adjusted)):
            if p_adj >= 0.05:
                continue
            
            group1, group2 = comp.split(' vs ')
            x1 = solidity_order.index(group1)
            x2 = solidity_order.index(group2)
            
            # Draw bar
            ax.plot([x1, x2], [bar_y, bar_y], 'k-', linewidth=1.5)
            ax.plot([x1, x1], [bar_y - bar_height/2, bar_y], 'k-', linewidth=1.5)
            ax.plot([x2, x2], [bar_y - bar_height/2, bar_y], 'k-', linewidth=1.5)
            
            # Add asterisks
            asterisks = get_asterisks(p_adj)
            ax.text((x1 + x2) / 2, bar_y + bar_height, asterisks, 
                   ha='center', va='bottom', fontsize=10, fontweight='bold')
            
            # Add p-value text
            ax.text((x1 + x2) / 2, bar_y + 2.5 * bar_height, f'p={p_adj:.3e}', 
                   ha='center', va='bottom', fontsize=8)
            
            bar_y += 3 * bar_height
        
        # Adjust ylim to accommodate annotations
        ax.set_ylim(bottom=prop_df['proportion_nonzero'].min() - 0.05 * y_range,
                   top=bar_y + 2 * bar_height)
    
    # Add grid
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_axisbelow(True)
    
    # Save figure
    os.makedirs(output_dir, exist_ok=True)
    output_pdf = os.path.join(output_dir, 'proportion_nonzero_boxplot.pdf')
    output_png = os.path.join(output_dir, 'proportion_nonzero_boxplot.png')
    
    plt.savefig(output_pdf, bbox_inches='tight', dpi=300)
    plt.savefig(output_png, bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"\n✅ Boxplot saved:")
    print(f"   - {output_pdf}")
    print(f"   - {output_png}")
    
    # Save statistics
    stats_df = pd.DataFrame({
        'Comparison': comparisons,
        'p_value': p_values,
        'p_adjusted': p_adjusted if len(p_values) > 0 else []
    })
    stats_path = os.path.join(output_dir, 'proportion_nonzero_stats.csv')
    stats_df.to_csv(stats_path, index=False)
    print(f"   - Statistics: {stats_path}")


# ============================================================================
# Main
# ============================================================================

def main():
    """Main function."""
    print("="*70)
    print("Spatial Multipage PDF Generation")
    print("="*70)
    
    # Load data
    print(f"\nLoading data from: {ADATA_PATH}")
    adata = sc.read(ADATA_PATH)
    print(f"   Loaded {len(adata)} spots, {adata.n_vars} genes")
    
    # Check for required columns
    print("\nChecking required columns...")
    required_cols = ['library_id', 'pathologist_annotation']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise KeyError(f"Missing required columns: {missing_cols}")
    print(f"   ✓ All required columns present")
    
    # Check for signature columns
    print("\nChecking signature columns...")
    if PROLIFERATION_SIG not in adata.obs.columns:
        print(f"   ⚠ Warning: {PROLIFERATION_SIG} not found")
        print(f"   Available signature columns (first 10):")
        sig_cols = [col for col in adata.obs.columns if col.startswith('sig:')]
        for col in sig_cols[:10]:
            print(f"      - {col}")
    else:
        print(f"   ✓ Found: {PROLIFERATION_SIG}")
    
    if MACROPHAGE_SIG not in adata.obs.columns:
        print(f"   ⚠ Warning: {MACROPHAGE_SIG} not found")
        print(f"   Available macrophage signatures:")
        mac_cols = [col for col in adata.obs.columns if 'macrophage' in col.lower() or 'mac' in col.lower()]
        for col in mac_cols[:10]:
            print(f"      - {col}")
    else:
        print(f"   ✓ Found: {MACROPHAGE_SIG}")
    
    # Prepare data for analysis
    print("\n" + "="*70)
    print("Preparing data...")
    print("="*70)
    
    # Prepare solidity
    adata.obs['solidity'] = prepare_solidity(adata)
    
    # Calculate truncated similarity
    print("Calculating truncated similarity...")
    if PROLIFERATION_SIG in adata.obs.columns and MACROPHAGE_SIG in adata.obs.columns:
        adata.obs['truncated_similarity'] = calculate_truncated_similarity(
            adata.obs[PROLIFERATION_SIG],
            adata.obs[MACROPHAGE_SIG]
        )
    else:
        print(f"Warning: Missing signature columns. Proliferation: {PROLIFERATION_SIG in adata.obs.columns}, "
              f"Macrophage: {MACROPHAGE_SIG in adata.obs.columns}")
        adata.obs['truncated_similarity'] = 0.0
    
    # Create multipage PDF
    print("\n" + "="*70)
    print("Creating multipage PDF...")
    print("="*70)
    create_multipage_pdf(adata, OUTPUT_PDF, PROLIFERATION_SIG, MACROPHAGE_SIG)
    
    # Create proportion boxplot
    plot_proportion_boxplot(adata, OUTPUT_DIR)
    
    print("\n" + "="*70)
    print("✅ Complete!")
    print("="*70)
    print(f"Output: {OUTPUT_PDF}")
    print(f"Boxplot: {os.path.join(OUTPUT_DIR, 'proportion_nonzero_boxplot.pdf')}")


if __name__ == "__main__":
    main()
