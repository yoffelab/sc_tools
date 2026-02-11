"""
Spatial Analysis: Spatially Variable Genes (SVG) and Ligand-Receptor Crosstalk.

1. SVG Discovery: Identify genes that characterize the transition from GGO to Solid tumors
   using Moran's I or squidpy.gr.spatial_autocorr.

2. Ligand-Receptor Analysis: Infer spatially constrained ligand-receptor interactions
   using squidpy.gr.ligrec.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import squidpy as sq
import json
import os
from pathlib import Path
from tqdm import tqdm

sc.settings.set_figure_params(dpi=300, dpi_save=400)

# Load data
print("Loading data...")
adata = sc.read('results/adata.normalized.scored.p35.h5ad')

# Ensure spatial coordinates are available
assert 'spatial' in adata.obsm, "Spatial coordinates not found in adata.obsm['spatial']"
print(f"Data shape: {adata.shape}")
print(f"Spatial coordinates shape: {adata.obsm['spatial'].shape}")

# Define tumor type mapping
pa = {
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

def keep_first_unique(input_list):
    seen = set()
    unique_list = []
    for item in input_list:
        if item not in seen:
            unique_list.append(item)
            seen.add(item)
    return unique_list

adata.obs['tumor_type'] = pd.Categorical(
    adata.obs['pathologist_annotation'].astype(str).replace(pa),
    categories=keep_first_unique(['Normal Alveolar Cells', 'Non-Solid Tumor', 'Solid Tumor'])
)

# Filter to tumor types for comparison
tumor_adata = adata[adata.obs['tumor_type'].isin(['Non-Solid Tumor', 'Solid Tumor'])].copy()
tumor_adata.obs['tumor_type'] = tumor_adata.obs['tumor_type'].cat.remove_unused_categories()

print(f"\nTumor type distribution:\n{tumor_adata.obs['tumor_type'].value_counts()}")

# ============================================================================
# 1. Spatially Variable Genes (SVG) Discovery
# ============================================================================
print("\n" + "="*70)
print("1. SPATIALLY VARIABLE GENES (SVG) DISCOVERY")
print("="*70)

os.makedirs('results/svg', exist_ok=True)
os.makedirs('figures/spatial/svg', exist_ok=True)

# Compute spatial neighbors graph
print("\nComputing spatial neighbor graph...")
sq.gr.spatial_neighbors(tumor_adata, coord_type="generic", delaunay=True)

# Compute spatial autocorrelation (Moran's I) for all genes
print("\nComputing spatial autocorrelation (Moran's I)...")
# Use a subset of highly variable genes to speed up computation
if 'highly_variable' not in tumor_adata.var.columns:
    print("Computing highly variable genes...")
    sc.pp.highly_variable_genes(tumor_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

hvg_genes = tumor_adata.var[tumor_adata.var['highly_variable']].index.tolist()
print(f"Computing Moran's I for {len(hvg_genes)} highly variable genes...")

sq.gr.spatial_autocorr(
    tumor_adata,
    mode="moran",
    n_perms=100,
    n_jobs=1,
    genes=hvg_genes[:500] if len(hvg_genes) > 500 else hvg_genes,  # Limit for speed
)

# Extract results
if 'moranI' in tumor_adata.uns:
    moran_results = tumor_adata.uns['moranI']
    moran_df = pd.DataFrame(moran_results).T
    moran_df.columns = ['I', 'pval_norm', 'pval_norm_fdr_bh']
    moran_df = moran_df.sort_values('I', ascending=False)
    
    # Save results
    moran_df.to_csv('results/svg/morans_i_results.csv')
    print(f"\nTop 20 spatially variable genes (by Moran's I):")
    print(moran_df.head(20))
    
    # Plot top SVGs
    top_svgs = moran_df.head(20).index.tolist()
    
    # Plot spatial distribution of top SVGs
    print("\nPlotting spatial distribution of top SVGs...")
    for lib_id in tumor_adata.obs['library_id'].unique()[:3]:  # Limit to first 3 libraries
        adata_lib = tumor_adata[tumor_adata.obs['library_id'] == lib_id].copy()
        
        for gene in top_svgs[:5]:  # Top 5 genes
            if gene in adata_lib.var_names:
                sc.pl.spatial(
                    adata_lib,
                    color=gene,
                    library_id=lib_id,
                    cmap='viridis',
                    show=False,
                    title=f'{gene} (Moran\'s I = {moran_df.loc[gene, "I"]:.3f})'
                )
                plt.savefig(
                    f'figures/spatial/svg/{lib_id}_{gene}_spatial.pdf',
                    bbox_inches='tight', dpi=300
                )
                plt.savefig(
                    f'figures/spatial/svg/{lib_id}_{gene}_spatial.png',
                    bbox_inches='tight', dpi=300
                )
                plt.close()
    
    # Compare SVG expression between Non-Solid and Solid
    print("\nComparing SVG expression between Non-Solid and Solid tumors...")
    top_10_svgs = moran_df.head(10).index.tolist()
    
    # Filter to genes present in data
    top_10_svgs = [g for g in top_10_svgs if g in tumor_adata.var_names]
    
    if len(top_10_svgs) > 0:
        # Create comparison plot
        fig, axes = plt.subplots(2, 5, figsize=(15, 6))
        axes = axes.flatten()
        
        for i, gene in enumerate(top_10_svgs[:10]):
            ax = axes[i]
            
            # Prepare data for boxplot
            plot_data = []
            labels = []
            for tumor_type in ['Non-Solid Tumor', 'Solid Tumor']:
                mask = tumor_adata.obs['tumor_type'] == tumor_type
                if mask.sum() > 0:
                    values = tumor_adata[mask, gene].X
                    if hasattr(values, 'toarray'):
                        values = values.toarray().flatten()
                    else:
                        values = values.flatten()
                    plot_data.append(values)
                    labels.append(tumor_type.replace(' Tumor', ''))
            
            if len(plot_data) == 2:
                bp = ax.boxplot(plot_data, labels=labels, patch_artist=True, showfliers=False)
                bp['boxes'][0].set_facecolor('#fc8d62')
                bp['boxes'][1].set_facecolor('#8da0cb')
                ax.set_title(f'{gene}\n(I={moran_df.loc[gene, "I"]:.3f})', fontsize=9)
                ax.set_ylabel('Expression', fontsize=8)
                ax.tick_params(labelsize=8)
        
        plt.suptitle('Top 10 Spatially Variable Genes: Non-Solid vs Solid', fontsize=12, fontweight='bold')
        plt.tight_layout()
        plt.savefig('figures/spatial/svg/top10_svg_comparison.pdf', bbox_inches='tight', dpi=300)
        plt.savefig('figures/spatial/svg/top10_svg_comparison.png', bbox_inches='tight', dpi=300)
        plt.close()

# ============================================================================
# 2. Ligand-Receptor Analysis
# ============================================================================
print("\n" + "="*70)
print("2. LIGAND-RECEPTOR CROSSTALK ANALYSIS")
print("="*70)

os.makedirs('results/ligrec', exist_ok=True)
os.makedirs('figures/spatial/ligrec', exist_ok=True)

# Ensure we have normalized data for ligand-receptor analysis
if tumor_adata.raw is not None:
    print("Using raw counts for normalization...")
    tumor_adata.X = tumor_adata.raw.X

# Normalize if needed (ligrec typically works on normalized data)
if np.max(tumor_adata.X) > 100:  # Likely raw counts
    print("Normalizing data for ligand-receptor analysis...")
    sc.pp.normalize_total(tumor_adata, target_sum=1e4)
    sc.pp.log1p(tumor_adata)

# Compute spatial neighbors if not already done
if 'spatial_neighbors' not in tumor_adata.obsp:
    print("Computing spatial neighbors for ligand-receptor analysis...")
    sq.gr.spatial_neighbors(tumor_adata, coord_type="generic", delaunay=True)

# Run ligand-receptor analysis
print("\nRunning ligand-receptor interaction analysis...")
print("This may take several minutes...")

try:
    sq.gr.ligrec(  # Note: may be sq.gr.ligprec in some squidpy versions
        tumor_adata,
        n_perms=1000,
        cluster_key='tumor_type',
        copy=True,
    )
    
    # Extract results
    if 'ligrec' in tumor_adata.uns:
        ligrec_results = tumor_adata.uns['ligrec']
        
        # Save results
        if isinstance(ligrec_results, dict) and 'pvalues' in ligrec_results:
            ligrec_df = ligrec_results['pvalues']
            ligrec_df.to_csv('results/ligrec/ligrec_pvalues.csv')
            print(f"\nLigand-receptor results saved to results/ligrec/ligrec_pvalues.csv")
            
            # Plot top interactions
            print("\nPlotting top ligand-receptor interactions...")
            
            # Get top interactions (lowest p-values)
            if hasattr(ligrec_df, 'stack'):
                top_interactions = ligrec_df.stack().nsmallest(20)
                print(f"\nTop 20 ligand-receptor interactions:")
                print(top_interactions)
            
            # Visualize interactions spatially for a sample library
            for lib_id in tumor_adata.obs['library_id'].unique()[:2]:  # First 2 libraries
                adata_lib = tumor_adata[tumor_adata.obs['library_id'] == lib_id].copy()
                
                # Plot spatial distribution of tumor types
                sc.pl.spatial(
                    adata_lib,
                    color='tumor_type',
                    library_id=lib_id,
                    show=False,
                    title=f'{lib_id}: Tumor Type Distribution'
                )
                plt.savefig(
                    f'figures/spatial/ligrec/{lib_id}_tumor_type_spatial.pdf',
                    bbox_inches='tight', dpi=300
                )
                plt.savefig(
                    f'figures/spatial/ligrec/{lib_id}_tumor_type_spatial.png',
                    bbox_inches='tight', dpi=300
                )
                plt.close()
        
        print("\n✅ Ligand-receptor analysis complete!")
        
except Exception as e:
    print(f"\n⚠️  Ligand-receptor analysis encountered an error: {e}")
    print("   This may be due to data format or missing dependencies.")
    print("   Continuing with SVG results...")

print("\n" + "="*70)
print("✅ SPATIAL ANALYSIS COMPLETE")
print("="*70)
print(f"   - SVG results: results/svg/morans_i_results.csv")
print(f"   - SVG plots: figures/spatial/svg/")
if os.path.exists('results/ligrec/ligrec_pvalues.csv'):
    print(f"   - Ligand-receptor results: results/ligrec/ligrec_pvalues.csv")
    print(f"   - Ligand-receptor plots: figures/spatial/ligrec/")

