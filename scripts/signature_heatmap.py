"""
Gene Signature Heatmap: Create heatmaps and clustermaps of signature scores.

This script creates heatmaps with:
- Rows: spots
- Columns: gene signatures
- Annotations: patient/sample (library_id) and solidity (tumor_type)

Produces 4 figures total:
1. Heatmap sorted by patient then solidity
2. Heatmap sorted by solidity then patient
3. Clustermap sorted by patient then solidity (with clustering applied after subsetting)
4. Clustermap sorted by solidity then patient (with clustering applied after subsetting)

For clustermaps, the script ensures large group boundaries are preserved while
allowing fine-scale clustering within groups.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
from pathlib import Path
import os
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

sc.settings.set_figure_params(dpi=300, dpi_save=400)

# ============================================================================
# Configuration
# ============================================================================

# Input/Output paths
ADATA_PATH = Path('results/adata.normalized.scored.p35.h5ad')
OUTPUT_DIR = Path('figures/manuscript/signature_heatmaps')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Column names
LIBRARY_ID_COL = 'library_id'
TUMOR_TYPE_COL = 'tumor_type'
SIGNATURE_PREFIX = 'sig:'
SIGNATURE_SUFFIX = '_z'

# Clustermap parameters
CLUSTER_METHOD = 'average'  # Linkage method
CLUSTER_METRIC = 'euclidean'  # Distance metric

# Heatmap parameters
HEATMAP_FIGSIZE = (16, 12)  # Will be adjusted based on data size
CLUSTERMAP_FIGSIZE = (18, 14)
DENDROGRAM_RATIO = 0.1


# ============================================================================
# Helper Functions
# ============================================================================

def get_signature_columns(adata, prefix=SIGNATURE_PREFIX, suffix=SIGNATURE_SUFFIX):
    """Get all signature columns from adata.obs."""
    sig_cols = [col for col in adata.obs.columns 
                if col.startswith(prefix) and col.endswith(suffix)]
    return sorted(sig_cols)


def prepare_tumor_type(adata):
    """Prepare tumor_type column if it doesn't exist."""
    if TUMOR_TYPE_COL not in adata.obs.columns:
        # Check if pathologist_annotation exists
        if 'pathologist_annotation' not in adata.obs.columns:
            raise KeyError(f"Neither '{TUMOR_TYPE_COL}' nor 'pathologist_annotation' found in adata.obs.columns")
        
        # Create tumor type mapping from pathologist_annotation
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
        
        adata.obs[TUMOR_TYPE_COL] = pd.Categorical(
            adata.obs['pathologist_annotation'].astype(str).replace(pa),
            categories=keep_first_unique(['Normal Alveolar Cells', 'Non-Solid Tumor', 'Solid Tumor'])
        )
    
    # Filter to only valid tumor types
    valid_types = ['Normal Alveolar Cells', 'Non-Solid Tumor', 'Solid Tumor']
    adata_sub = adata[adata.obs[TUMOR_TYPE_COL].isin(valid_types)].copy()
    
    # Only remove unused categories if it's a categorical
    if isinstance(adata_sub.obs[TUMOR_TYPE_COL].dtype, pd.CategoricalDtype):
        adata_sub.obs[TUMOR_TYPE_COL] = adata_sub.obs[TUMOR_TYPE_COL].cat.remove_unused_categories()
    
    return adata_sub


def create_annotation_dataframe(adata_sub):
    """Create annotation dataframe for heatmap."""
    annotations = pd.DataFrame({
        'Patient': adata_sub.obs[LIBRARY_ID_COL].values,
        'Solidity': adata_sub.obs[TUMOR_TYPE_COL].values
    }, index=adata_sub.obs_names)
    
    # Create shortened labels for display
    annotations['Patient_short'] = annotations['Patient'].astype(str)
    annotations['Solidity_short'] = annotations['Solidity'].str.replace(' Alveolar Cells', '').str.replace(' Tumor', '')
    
    return annotations


def hex_to_rgb(hex_color):
    """Convert hex color to RGB tuple (0-1 range)."""
    hex_color = hex_color.lstrip('#')
    if len(hex_color) != 6:
        return (0.8, 0.8, 0.8)  # Default gray
    return tuple(int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4))


def create_colors_dict(annotations):
    """Create color mapping dictionaries for annotations.
    
    Returns colors as RGB tuples in (0-1) range for matplotlib.
    """
    # Patient colors (using qualitative palette) - already RGB tuples
    unique_patients = sorted(annotations['Patient'].unique())
    patient_colors_rgb = sns.color_palette("Set3", len(unique_patients))
    patient_colors = dict(zip(unique_patients, patient_colors_rgb))
    
    # Solidity colors (default mapping as hex, convert to RGB)
    default_solidity_colors_hex = {
        'Normal Alveolar Cells': '#66c2a5',
        'Non-Solid Tumor': '#fc8d62',
        'Solid Tumor': '#8da0cb',
    }
    
    # Get unique solidity values and assign colors
    unique_solidities = sorted(annotations['Solidity'].unique())
    solidity_colors = {}
    for solidity in unique_solidities:
        if solidity in default_solidity_colors_hex:
            solidity_colors[solidity] = hex_to_rgb(default_solidity_colors_hex[solidity])
        else:
            # Assign a gray color for unknown solidity types
            solidity_colors[solidity] = (0.8, 0.8, 0.8)  # RGB gray
            print(f"   Warning: Unknown solidity type '{solidity}', assigning gray color")
    
    return patient_colors, solidity_colors


def cluster_within_groups(data_matrix, group_labels, method=CLUSTER_METHOD, metric=CLUSTER_METRIC):
    """
    Cluster data while preserving group boundaries.
    
    Parameters
    ----------
    data_matrix : np.ndarray
        Data matrix (spots x signatures)
    group_labels : np.ndarray
        Group labels for each spot
    method : str
        Linkage method
    metric : str
        Distance metric
    
    Returns
    -------
    np.ndarray
        Reordered indices that cluster within groups but preserve group order
    """
    unique_groups = np.unique(group_labels)
    reordered_indices = []
    
    for group in unique_groups:
        group_mask = group_labels == group
        group_indices = np.where(group_mask)[0]
        group_data = data_matrix[group_indices, :]
        
        if len(group_indices) > 1:
            # Cluster within this group
            distances = pdist(group_data, metric=metric)
            if len(distances) > 0:
                linkage_matrix = linkage(distances, method=method)
                leaves = leaves_list(linkage_matrix)
                # Map back to original indices
                group_clustered = group_indices[leaves]
            else:
                group_clustered = group_indices
        else:
            group_clustered = group_indices
        
        reordered_indices.extend(group_clustered.tolist())
    
    return np.array(reordered_indices)


# ============================================================================
# Main Functions
# ============================================================================

def create_signature_heatmap(adata_sub, sig_columns, sort_by=['patient', 'solidity'], 
                            output_suffix='patient_solidity', cluster=False):
    """
    Create heatmap or clustermap of signature scores.
    
    Parameters
    ----------
    adata_sub : AnnData
        Annotated data object with signature scores
    sig_columns : list
        List of signature column names
    sort_by : list
        Sorting order: ['patient', 'solidity'] or ['solidity', 'patient']
    output_suffix : str
        Suffix for output filename
    cluster : bool
        If True, create clustermap. If False, create regular heatmap.
    """
    # Extract signature scores
    sig_matrix = adata_sub.obs[sig_columns].values
    sig_df = pd.DataFrame(sig_matrix, 
                         index=adata_sub.obs_names, 
                         columns=[s.replace(SIGNATURE_PREFIX, '').replace(SIGNATURE_SUFFIX, '') 
                                 for s in sig_columns])
    
    # Create annotations
    annotations = create_annotation_dataframe(adata_sub)
    
    # Check that annotations were created correctly
    if 'Patient' not in annotations.columns or 'Solidity' not in annotations.columns:
        raise KeyError(f"Annotations dataframe missing required columns. Found: {annotations.columns.tolist()}")
    
    patient_colors, solidity_colors = create_colors_dict(annotations)
    
    # Define solidity order (always: Normal, Non-Solid, Solid)
    # Make it accessible for legend creation
    solidity_order = ['Normal Alveolar Cells', 'Non-Solid Tumor', 'Solid Tumor']
    
    # Calculate mean scores per patient for ordering
    patient_means = {}
    unique_patients = annotations['Patient'].unique()
    for patient in unique_patients:
        patient_mask = annotations['Patient'] == patient
        patient_indices = annotations.index[patient_mask]
        patient_sig_data = sig_df.loc[patient_indices]
        patient_means[patient] = patient_sig_data.mean(axis=1).mean()
    
    # Order patients by mean scores (highest first)
    patient_order = sorted(unique_patients, key=lambda x: patient_means[x], reverse=True)
    
    # Create grouping keys based on sort_by
    if sort_by == ['patient', 'solidity']:
        # Group by: (patient, solidity)
        annotations['group_key'] = annotations.apply(
            lambda row: (patient_order.index(row['Patient']), 
                        solidity_order.index(row['Solidity'])), 
            axis=1
        )
    elif sort_by == ['solidity', 'patient']:
        # Group by: (solidity, patient)
        annotations['group_key'] = annotations.apply(
            lambda row: (solidity_order.index(row['Solidity']),
                        patient_order.index(row['Patient'])), 
            axis=1
        )
    else:
        raise ValueError(f"Invalid sort_by: {sort_by}")
    
    # Sort by group key first (this establishes the macro order)
    annotations_sorted = annotations.sort_values('group_key')
    sorted_indices = annotations_sorted.index
    
    # Reorder data
    sig_df_sorted = sig_df.loc[sorted_indices]
    annotations_sorted = annotations_sorted.drop('group_key', axis=1)
    
    if cluster:
        # For clustermap: Split by the 2 layers, cluster within each group
        if sort_by == ['patient', 'solidity']:
            group_cols = ['Patient', 'Solidity']
        else:  # ['solidity', 'patient']
            group_cols = ['Solidity', 'Patient']
        
        # Get unique groups (observed=True to silence warning)
        grouped = annotations_sorted.groupby(group_cols, observed=True)
        clustered_indices_list = []
        
        for group_key, group_df in grouped:
            group_indices = group_df.index.tolist()
            
            if len(group_indices) > 1:
                # Cluster within this group
                group_data = sig_df_sorted.loc[group_indices].values
                # cluster_within_groups returns positions within the group
                # All items are in the same group, so just pass a dummy label
                clustered_positions = cluster_within_groups(
                    group_data,
                    np.zeros(len(group_indices), dtype=int),  # All same group (dummy)
                    method=CLUSTER_METHOD,
                    metric=CLUSTER_METRIC
                )
                # Map positions back to original indices
                clustered_indices_list.extend([group_indices[pos] for pos in clustered_positions])
            else:
                clustered_indices_list.extend(group_indices)
        
        # Apply clustering
        sig_df_clustered = sig_df_sorted.loc[clustered_indices_list]
        annotations_clustered = annotations_sorted.loc[clustered_indices_list]
    else:
        sig_df_clustered = sig_df_sorted
        annotations_clustered = annotations_sorted
    
    # Prepare annotation colors as RGB tuples
    # Handle missing keys gracefully with gray RGB tuple
    default_gray = (0.8, 0.8, 0.8)
    patient_color_list = [patient_colors.get(p, default_gray) for p in annotations_clustered['Patient']]
    solidity_color_list = [solidity_colors.get(s, default_gray) for s in annotations_clustered['Solidity']]
    
    # Check for missing colors
    missing_patients = set(annotations_clustered['Patient']) - set(patient_colors.keys())
    missing_solidities = set(annotations_clustered['Solidity']) - set(solidity_colors.keys())
    if missing_patients:
        print(f"   Warning: Missing colors for patients: {missing_patients}")
    if missing_solidities:
        print(f"   Warning: Missing colors for solidities: {missing_solidities}")
    
    # Convert to numpy array with proper shape for imshow: (1, n_spots, 3)
    patient_colors_array = np.array(patient_color_list).reshape(1, -1, 3)
    solidity_colors_array = np.array(solidity_color_list).reshape(1, -1, 3)
    
    annotation_colors = pd.DataFrame({
        'Patient': patient_color_list,
        'Solidity': solidity_color_list
    }, index=annotations_clustered.index)
    
    if cluster:
        # Create clustermap
        g = sns.clustermap(
            sig_df_clustered.T,  # Transpose: signatures as rows, spots as columns
            cmap='RdBu_r',
            center=0,
            vmin=-3,
            vmax=3,
            figsize=CLUSTERMAP_FIGSIZE,
            linewidths=0,
            cbar_kws={'label': 'Z-scored Signature Score'},
            dendrogram_ratio=DENDROGRAM_RATIO,
            cbar_pos=(0.02, 0.85, 0.03, 0.12),
            row_cluster=True,
            col_cluster=False,  # Don't cluster columns (already sorted/clustered within groups)
            yticklabels=True,
            xticklabels=False,  # Too many spots to show labels
            method=CLUSTER_METHOD,
            metric=CLUSTER_METRIC
        )
        
        # Add column annotations
        # Note: clustermap transposes, so we need to add annotations to columns
        # But our annotations are for spots (rows in original data, columns after transpose)
        # The clustermap already handles this correctly
        
        # Set title
        title = f'Signature Scores Clustermap (sorted by {sort_by[0]} then {sort_by[1]})'
        g.fig.suptitle(title, fontsize=14, fontweight='bold', y=0.995)
        
        # Add annotation colors manually (clustermap doesn't support col_colors directly)
        # Get heatmap position
        heatmap_pos = g.ax_heatmap.get_position()
        
        # Create annotation bars
        n_spots = len(annotations_clustered)
        
        # Calculate annotation bar height
        ann_height = 0.015
        
        # Add patient annotation bar (above heatmap)
        ax_ann_patient = g.fig.add_axes([
            heatmap_pos.x0,
            heatmap_pos.y1,
            heatmap_pos.width,
            ann_height
        ])
        # patient_colors_array was created earlier with shape (1, n_spots, 3)
        ax_ann_patient.imshow(patient_colors_array, aspect='auto', interpolation='nearest', rasterized=True)
        ax_ann_patient.set_xticks([])
        ax_ann_patient.set_yticks([])
        ax_ann_patient.set_ylabel('Patient', rotation=0, ha='right', va='center', fontsize=9)
        ax_ann_patient.spines['top'].set_visible(False)
        ax_ann_patient.spines['right'].set_visible(False)
        ax_ann_patient.spines['bottom'].set_visible(False)
        ax_ann_patient.spines['left'].set_visible(False)
        
        # Add solidity annotation bar (above patient)
        ax_ann_solidity = g.fig.add_axes([
            heatmap_pos.x0,
            heatmap_pos.y1 + ann_height,
            heatmap_pos.width,
            ann_height
        ])
        # solidity_colors_array was created earlier with shape (1, n_spots, 3)
        ax_ann_solidity.imshow(solidity_colors_array, aspect='auto', interpolation='nearest', rasterized=True)
        ax_ann_solidity.set_xticks([])
        ax_ann_solidity.set_yticks([])
        ax_ann_solidity.set_ylabel('Solidity', rotation=0, ha='right', va='center', fontsize=9)
        ax_ann_solidity.spines['top'].set_visible(False)
        ax_ann_solidity.spines['right'].set_visible(False)
        ax_ann_solidity.spines['bottom'].set_visible(False)
        ax_ann_solidity.spines['left'].set_visible(False)
        
        # Adjust layout
        g.ax_heatmap.set_xlabel(f'Spots (n={n_spots})', fontsize=11)
        g.ax_heatmap.set_ylabel('Gene Signatures', fontsize=11)
        
        # Add solidity legend in top right corner
        legend_elements = [
            Patch(facecolor=solidity_colors.get(s, (0.8, 0.8, 0.8)), 
                  label=s.replace(' Alveolar Cells', '').replace(' Tumor', ''))
            for s in solidity_order if s in solidity_colors
        ]
        legend_ax = g.fig.add_axes([0.95, 0.95, 0.05, 0.05])
        legend_ax.axis('off')
        legend = legend_ax.legend(handles=legend_elements, loc='upper right', 
                                 frameon=True, fontsize=9, title='Solidity', 
                                 title_fontsize=10)
        legend.get_frame().set_facecolor('white')
        legend.get_frame().set_alpha(0.9)
        
        plt.savefig(OUTPUT_DIR / f'signature_clustermap_{output_suffix}.pdf', 
                   bbox_inches='tight', dpi=300)
        plt.savefig(OUTPUT_DIR / f'signature_clustermap_{output_suffix}.png', 
                   bbox_inches='tight', dpi=300)
        plt.close()
        
    else:
        # Create regular heatmap
        fig, ax = plt.subplots(figsize=HEATMAP_FIGSIZE)
        
        sns.heatmap(
            sig_df_clustered.T,  # Transpose: signatures as rows, spots as columns
            cmap='RdBu_r',
            center=0,
            vmin=-3,
            vmax=3,
            cbar_kws={'label': 'Z-scored Signature Score'},
            linewidths=0,
            ax=ax,
            xticklabels=False,
            yticklabels=True
        )
        
        # Add annotations
        # Create annotation bars above heatmap
        n_spots = len(annotations_clustered)
        
        # Patient annotation
        ax_ann_patient = fig.add_axes([ax.get_position().x0,
                                      ax.get_position().y1,
                                      ax.get_position().width,
                                      0.02])
        ax_ann_patient.imshow(patient_colors_array, aspect='auto', interpolation='nearest')
        ax_ann_patient.set_xticks([])
        ax_ann_patient.set_yticks([])
        ax_ann_patient.set_ylabel('Patient', rotation=0, ha='right', va='center', fontsize=10)
        
        # Solidity annotation
        ax_ann_solidity = fig.add_axes([ax.get_position().x0,
                                       ax.get_position().y1 + 0.02,
                                       ax.get_position().width,
                                       0.02])
        ax_ann_solidity.imshow(solidity_colors_array, aspect='auto', interpolation='nearest')
        ax_ann_solidity.set_xticks([])
        ax_ann_solidity.set_yticks([])
        ax_ann_solidity.set_ylabel('Solidity', rotation=0, ha='right', va='center', fontsize=10)
        
        ax.set_title(f'Signature Scores Heatmap (sorted by {sort_by[0]} then {sort_by[1]})', 
                    fontsize=14, fontweight='bold', pad=20)
        ax.set_xlabel(f'Spots (n={n_spots})', fontsize=12)
        ax.set_ylabel('Gene Signatures', fontsize=12)
        
        # Add solidity legend in top right corner
        legend_elements = [
            Patch(facecolor=solidity_colors.get(s, (0.8, 0.8, 0.8)), 
                  label=s.replace(' Alveolar Cells', '').replace(' Tumor', ''))
            for s in solidity_order if s in solidity_colors
        ]
        legend_ax = fig.add_axes([0.95, 0.95, 0.05, 0.05])
        legend_ax.axis('off')
        legend = legend_ax.legend(handles=legend_elements, loc='upper right', 
                                 frameon=True, fontsize=9, title='Solidity', 
                                 title_fontsize=10)
        legend.get_frame().set_facecolor('white')
        legend.get_frame().set_alpha(0.9)
        
        plt.savefig(OUTPUT_DIR / f'signature_heatmap_{output_suffix}.pdf', 
                   bbox_inches='tight', dpi=300)
        plt.savefig(OUTPUT_DIR / f'signature_heatmap_{output_suffix}.png', 
                   bbox_inches='tight', dpi=300)
        plt.close()
    
    print(f"   ✅ Saved: signature_{'clustermap' if cluster else 'heatmap'}_{output_suffix}.pdf")


# ============================================================================
# Main Workflow
# ============================================================================

def main():
    """Main function to generate signature heatmaps."""
    print("="*70)
    print("GENE SIGNATURE HEATMAP GENERATION")
    print("="*70)
    
    # Load data
    print("\n1. Loading data...")
    if not ADATA_PATH.exists():
        raise FileNotFoundError(f"AnnData file not found: {ADATA_PATH}")
    
    adata = sc.read_h5ad(ADATA_PATH)
    print(f"   Loaded: {adata.shape} (spots x genes)")
    
    # Prepare tumor type
    print("\n2. Preparing tumor type annotations...")
    adata_sub = prepare_tumor_type(adata)
    print(f"   Filtered to {adata_sub.n_obs} spots with valid tumor types")
    print(f"   Tumor type distribution:")
    print(adata_sub.obs[TUMOR_TYPE_COL].value_counts())
    
    # Get signature columns
    print("\n3. Extracting signature columns...")
    sig_columns = get_signature_columns(adata_sub)
    print(f"   Found {len(sig_columns)} signature columns")
    
    if len(sig_columns) == 0:
        print("   ⚠ No signature columns found. Please run score_gene_signatures.py first.")
        return
    
    # Check for library_id
    if LIBRARY_ID_COL not in adata_sub.obs.columns:
        print(f"   ⚠ {LIBRARY_ID_COL} column not found. Cannot create patient annotations.")
        return
    
    # Generate heatmaps
    print("\n4. Generating heatmaps...")
    
    # 1. Heatmap sorted by patient then solidity
    print("\n   Creating heatmap (patient → solidity)...")
    create_signature_heatmap(adata_sub, sig_columns, 
                            sort_by=['patient', 'solidity'],
                            output_suffix='patient_solidity',
                            cluster=False)
    
    # 2. Heatmap sorted by solidity then patient
    print("\n   Creating heatmap (solidity → patient)...")
    create_signature_heatmap(adata_sub, sig_columns,
                            sort_by=['solidity', 'patient'],
                            output_suffix='solidity_patient',
                            cluster=False)
    
    # 3. Clustermap sorted by patient then solidity
    print("\n   Creating clustermap (patient → solidity)...")
    create_signature_heatmap(adata_sub, sig_columns,
                            sort_by=['patient', 'solidity'],
                            output_suffix='patient_solidity',
                            cluster=True)
    
    # 4. Clustermap sorted by solidity then patient
    print("\n   Creating clustermap (solidity → patient)...")
    create_signature_heatmap(adata_sub, sig_columns,
                            sort_by=['solidity', 'patient'],
                            output_suffix='solidity_patient',
                            cluster=True)
    
    print("\n" + "="*70)
    print("✅ SIGNATURE HEATMAP GENERATION COMPLETE")
    print("="*70)
    print(f"   All figures saved to: {OUTPUT_DIR}")
    print("\n   Generated files:")
    print("     - signature_heatmap_patient_solidity.pdf")
    print("     - signature_heatmap_solidity_patient.pdf")
    print("     - signature_clustermap_patient_solidity.pdf")
    print("     - signature_clustermap_solidity_patient.pdf")
    print("="*70)


if __name__ == "__main__":
    main()
