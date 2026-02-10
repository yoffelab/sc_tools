import scanpy as sc
import matplotlib.pyplot as plt
import sthd
import pandas as pd
import os
import squidpy as sq
import numpy as np
import warnings
import scipy.sparse as sp
from matplotlib.backends.backend_pdf import PdfPages

# Create figure directory only once
os.makedirs('figures', exist_ok=True)


def compute_and_plot_markers(
    adata,
    groupby_key,
    min_cells=1000,
    rank_method='wilcoxon',
    top_n_for_dict=5,
    n_genes_for_plots=10,
    output_dir='figures'
):
    """
    Full pipeline to compute markers, subset large clusters, and generate multiple plots.
    Returns: marker_dict {group: list of top N genes}
    """
    
    # 1. Compute marker genes
    sc.tl.rank_genes_groups(adata, groupby_key, method=rank_method, use_raw=False)
    
    # 2. Extract marker dataframe
    markers_df = sc.get.rank_genes_groups_df(adata, group=None)
    
    # 3. Print quick view of top 3 markers per group
    print(f"\nTop markers for {groupby_key}:")
    print(markers_df.groupby('group').head(3))
    
    # 4. Subset for larger clusters
    valid_groups = (
        adata.obs[groupby_key]
        .value_counts()
        .loc[lambda x: x > min_cells]
        .index.tolist()
    )
    adata_subset = adata[adata.obs[groupby_key].isin(valid_groups)].copy()
    adata_subset.obs[groupby_key] = pd.Categorical(adata_subset.obs[groupby_key], categories=valid_groups)
    
    # 5. Basic rank_genes_groups plot
    sc.pl.rank_genes_groups(adata, n_genes=n_genes_for_plots, sharey=False, show=False)
    plt.savefig(os.path.join(output_dir, f'{groupby_key}_signature_genes.pdf'))
    plt.close()

    # 6. Heatmap for top N genes per group
    sc.pl.rank_genes_groups_heatmap(
        adata_subset,
        groupby=groupby_key,
        n_genes=n_genes_for_plots,
        use_raw=False,
        show_gene_labels=True,
        swap_axes=True,
        show=False
    )
    plt.savefig(os.path.join(output_dir, f'genes_{groupby_key}_heatmap.pdf'))
    plt.close()

    # 7. Build top marker dictionary for future plotting
    top_genes_df = (
        markers_df.groupby("group")
        .apply(lambda x: x.nlargest(top_n_for_dict, "logfoldchanges"))
        .reset_index(drop=True)
    )
    marker_dict = (
        top_genes_df.groupby("group")['names']
        .apply(list)
        .to_dict()
    )

    # 8. Generate all downstream plots
    for plot_func, suffix in [
        (sc.pl.heatmap, "heatmap"),
        (sc.pl.tracksplot, "tracksplot"),
        (sc.pl.dotplot, "dotplot"),
        (sc.pl.matrixplot, "matrixplot"),
    ]:
        kwargs = {
            'var_names': marker_dict,
            'groupby': groupby_key,
            'use_raw': False,
            'show': False,
        }
        
        # Add arguments depending on plot type
        if suffix == "heatmap":
            kwargs.update({
                'swap_axes': True,
                'show_gene_labels': True,
                'standard_scale': 'var'
            })
        elif suffix == "matrixplot":
            kwargs.update({
                'swap_axes': True,
                'standard_scale': 'var'
            })
        elif suffix == "dotplot":
            kwargs.update({
                'standard_scale': 'var'
            })
        elif suffix == "tracksplot":
            kwargs.update({
                'dendrogram': False
            })

        plot_func(adata, **kwargs)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'top_genes_{groupby_key}_{suffix}.pdf'))
        plt.close()


    return marker_dict

adata = sc.read('output/scvi.leiden.grouped.h5ad')

for col in ['meta_cluster', 'sample', 'tissue_type']:
    sthd.pl.umap(adata, color = col, save_name = f'umap_{col}.png')

# adata = adata[:,~adata.var.index.str.contains('MT')]



# ================================
# 1. CRC gene sets
# ================================
gene_sets = {
    'epithelial': [
        'EPCAM', 'KRT8', 'KRT18', 'KRT20', 'CDH1', 'MUC2', 'MUC13', 
        'VIL1', 'FABP1', 'TFF3', 'REG4', 'CA1', 'CA2', 'CA4',
        'CEACAM5', 'LGR5', 'OLFM4', 'CDX2', 'SPINK1', 'SI'
    ],
    'tumor': [
        'MKI67', 'PCNA', 'TOP2A', 'BIRC5', 'MYC', 'AXIN2', 'CTNNB1',
        'CCND1', 'CDKN2A', 'CDKN1A', 'KRAS', 'TP53', 'PIK3CA',
        'SMAD4', 'FBXW7', 'BRAF', 'APC', 'PTEN', 'EGFR', 'ERBB2',
        'MCM2', 'MCM5', 'MCM6', 'CENPF', 'NEK2', 'CDCA3'
    ],
    'stromal': [
        'COL1A1', 'COL1A2', 'COL3A1', 'COL5A1', 'COL6A3', 'FN1',
        'ACTA2', 'TAGLN', 'PDGFRA', 'PDGFRB', 'POSTN', 'THBS1', 
        'SPARC', 'LOX', 'SERPINE1', 'TIMP1', 'MMP2', 'MMP9', 'MMP14',
        'FAP', 'CTHRC1', 'LUM', 'DCN', 'VCAN', 'VIM', 'S100A4'
    ],
    'immune': [
        'PTPRC', 'CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B',
        'FOXP3', 'PDCD1', 'CTLA4', 'CD79A', 'CD79B', 'MS4A1',
        'IGHM', 'IGKC', 'IGHA1',
        'CD68', 'CSF1R', 'CD163', 'MRC1', 'ITGAM',
        'S100A8', 'S100A9', 'LYZ', 'CXCL9', 'CXCL10', 'CCL5', 'IFNG',
        'NKG7', 'GZMB', 'PRF1', 'GNLY', 'FCGR3A'
    ]
}

# ================================
# 2. Filter gene sets based on your data
# ================================

# adata.var_names = adata.var_names.str.upper()
# available_genes = set(adata.var_names)

# for key in gene_sets:
#     gene_sets[key] = [gene for gene in gene_sets[key] if gene in available_genes]
#     print(f'{key}: {len(gene_sets[key])} genes retained')
sthd.pp.gene_set_score(adata, gene_sets)

# ================================
# 3. Calculate gene scores
# ================================

# for label, genes in gene_sets.items():
#     if len(genes) == 0:
#         warnings.warn(f"Skipping {label} — no genes found.")
#         continue
#     sc.tl.score_genes(adata, gene_list=genes, score_name=f'{label}_score', use_raw=False)

# ================================
# 4. Apply spatial smoothing per sample
# ================================
sthd.pp.gene_set_smoothen(adata, gene_sets)

# for label in gene_sets.keys():
#     colname = f'{label}_score_smoothed'
#     if colname not in adata.obs.columns:
#         adata.obs[colname] = np.nan

# for sample in adata.obs['sample'].unique():
#     print(f'\nProcessing sample: {sample}')
#     ad_sample = adata[adata.obs['sample'] == sample].copy()

#     sq.gr.spatial_neighbors(ad_sample, coord_type='generic')
#     A = ad_sample.obsp['spatial_connectivities']

#     for label in gene_sets.keys():
#         score_col = f'{label}_score'
#         if score_col in ad_sample.obs.columns:
#             scores = ad_sample.obs[score_col].values
#             smoothed = A.dot(scores)
#             smoothed = np.array(smoothed).flatten()
#             norm = np.array(A.sum(axis=1)).flatten()
#             smoothed = smoothed / np.where(norm == 0, 1, norm)

#             # Fully safe positional assignment
#             mask = adata.obs['sample'] == sample
#             positions = np.where(mask)[0]
#             adata.obs.iloc[positions, adata.obs.columns.get_loc(f'{score_col}_smoothed')] = smoothed

# ================================
# 5. Plotting (2x4 layout, fully stable)
# ================================
# Open one multipage PDF
with PdfPages('figures/all_samples_scores_panel.pdf') as pdf:

    for sample in adata.obs['sample'].unique():
        ad_sample = adata[adata.obs['sample'] == sample].copy()

        fig, axs = plt.subplots(2, 4, figsize=(20, 10), dpi = 150)

        for j, label in enumerate(gene_sets.keys()):
            score_col = f'{label}_score'
            score_col_smooth = f'{label}_score_smoothed'

            vmin = ad_sample.obs[[score_col, score_col_smooth]].min().min()
            vmax = ad_sample.obs[[score_col, score_col_smooth]].max().max()

            # Raw
            sc.pl.spatial(
                ad_sample,
                color=score_col,
                spot_size=1.5,
                cmap='viridis',
                vmin=vmin, vmax=vmax,
                show=False, ax=axs[0, j],
                title=f"{label} (raw)"
            )

            # Smoothed
            sc.pl.spatial(
                ad_sample,
                color=score_col_smooth,
                spot_size=1.5,
                cmap='viridis',
                vmin=vmin, vmax=vmax,
                show=False, ax=axs[1, j],
                title=f"{label} (smoothed)"
            )

        plt.suptitle(f'Sample: {sample}', fontsize=18)
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        pdf.savefig(fig)
        plt.close(fig)


# ================================
# 6. Softmax Labeling of Entity
# ================================
import numpy as np

# 1. Collect all score columns (you can include smoothed or raw depending on your choice)
score_cols = [col for col in adata.obs.columns if col.endswith('_score_smoothed')]

# 2. Extract raw scores as matrix
scores = adata.obs[score_cols].values

# 3. Normalize (z-score across spots for each score column)
means = scores.mean(axis=0, keepdims=True)
stds = scores.std(axis=0, keepdims=True)
normalized_scores = (scores - means) / np.where(stds == 0, 1, stds)

# 4. Compute softmax
def softmax(x):
    e_x = np.exp(x - np.max(x, axis=1, keepdims=True))  # numerical stability
    return e_x / e_x.sum(axis=1, keepdims=True)

softmax_scores = softmax(normalized_scores)

# 5. Store softmax scores as new columns
for i, label in enumerate(score_cols):
    adata.obs[f"{label}_softmax"] = softmax_scores[:, i]

# 6. Assign entity based on argmax
entity_labels = [label.replace('_score', '') for label in score_cols]
adata.obs['softmax_entity'] = softmax_scores.argmax(axis=1)
adata.obs['softmax_entity_label'] = adata.obs['softmax_entity'].map(lambda x: entity_labels[x])

# ================================
# 7. Softmax Labeling Plotting
# ================================
from matplotlib.backends.backend_pdf import PdfPages

with PdfPages('figures/all_samples_softmax_panel.pdf') as pdf:
    for sample in adata.obs['sample'].unique():
        ad_sample = adata[adata.obs['sample'] == sample].copy()
        fig, axs = plt.subplots(2, 3, figsize=(15, 10))

        softmax_cols = [col for col in ad_sample.obs.columns if col.endswith('_softmax')]
        entity_labels = [col.replace('_score_softmax', '') for col in softmax_cols]

        # Plot softmax scores on axs[0, 0], axs[0, 1], axs[0, 2], axs[1, 0]
        softmax_positions = [(0,0), (0,1), (0,2), (1,0)]
        for (softmax_col, entity_label), (r, c) in zip(zip(softmax_cols, entity_labels), softmax_positions):
            sc.pl.spatial(
                ad_sample,
                color=softmax_col,
                spot_size=1.5,
                cmap='viridis',
                vmin=0, vmax=1,
                show=False, ax=axs[r, c],
                title=f"{entity_label} (softmax)"
            )

        # Assigned entity label on axs[1,1]
        sc.pl.spatial(
            ad_sample,
            color="softmax_entity_label",
            spot_size=1.5,
            cmap='tab20',
            show=False, ax=axs[1, 1],
            title="Softmax Assigned Label"
        )

        # Disable empty last panel axs[1,2]
        axs[1, 2].axis('off')

        plt.suptitle(f'Sample: {sample}', fontsize=18)
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        pdf.savefig(fig)
        plt.close(fig)





# --- Apply function to multiple groupings ---
groupings = ["tissue_type", "meta_cluster", "softmax_entity_label"]
all_marker_dicts = {}

for group in groupings:
    marker_dict = compute_and_plot_markers(
        adata=adata,
        groupby_key=group,
        min_cells=1000,
        rank_method='wilcoxon',
        top_n_for_dict=5,
        n_genes_for_plots=10,
        output_dir='figures'
    )
    all_marker_dicts[group] = marker_dict


fig, ax = plt.subplots(1,1,figsize = (3,3), dpi = 300)
sc.pl.umap(adata, color = ['softmax_entity_label'], save = 'umap_softmax_entity_label.png', show = False, ax = ax)


import squidpy as sq
# Compute spatial neighbors
sq.gr.spatial_neighbors(adata)
sq.gr.spatial_autocorr(adata, mode='moran')
sq.pl.spatial_autocorr(adata, mode='moran', genes=adata.var_names[:4])
