import scanpy as sc
import tangram as tg
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from dask import delayed, compute
from dask.diagnostics import ProgressBar

# -------------------------------
# Configuration
# -------------------------------
sc_data_file = "results/seurat_object.h5ad"
visium_hd_data_file = "results/adata.annotation.masked.h5ad"
celltype_key = 'cell.type'
visium_batch_key = 'sample_id'
sc_batch_key = 'Batch'
output_dir = 'output/tangram_batches'
os.makedirs(output_dir, exist_ok=True)
n_genes_max = 2000

# -------------------------------
# Load data
# -------------------------------
sc_adata = sc.read_h5ad(sc_data_file)
sc_adata.var_names_make_unique()

visium_adata = sc.read_h5ad(visium_hd_data_file)
visium_adata.var_names_make_unique()
visium_adata.obs_names_make_unique()

# -------------------------------
# scRNA-seq: log-normalize if needed
# -------------------------------
if np.max(sc_adata.X) > 100:  # crude check for raw counts
    print("Normalizing and log-transforming scRNA-seq data...")
    sc.pp.normalize_total(sc_adata, target_sum=1e4)
    sc.pp.log1p(sc_adata)
else:
    print("scRNA-seq data already appears log-normalized")

# -------------------------------
# Visium: use raw counts if available
# -------------------------------
if visium_adata.raw is not None:
    print("Using visium raw counts")
    raw_adata = visium_adata.raw.to_adata()
    visium_adata.X = raw_adata[:, visium_adata.var_names].X

# Check if Visium data is count-like
if hasattr(visium_adata.X, "data"):
    assert np.allclose(visium_adata.X.data.astype(int), visium_adata.X.data), "visium_adata.X contains non-integers"
else:
    assert np.allclose(visium_adata.X.astype(int), visium_adata.X), "visium_adata.X contains non-integers"

# -------------------------------
# Signature gene selection (~2000 genes)
# -------------------------------
print("Selecting signature genes...")

# Step 1: HVGs
sc.pp.highly_variable_genes(sc_adata, flavor="seurat_v3", n_top_genes=3000, subset=False, batch_key=sc_batch_key)
hvg_genes = set(sc_adata.var[sc_adata.var.highly_variable].index)

# Step 2: Top markers per cell type
sc.tl.rank_genes_groups(sc_adata, groupby=celltype_key, method='wilcoxon', use_raw=False)
marker_genes = []
for group in sc_adata.obs[celltype_key].unique():
    top_genes = sc.get.rank_genes_groups_df(sc_adata, group=group).head(100).names.tolist()
    marker_genes.extend(top_genes)
marker_genes = set(marker_genes)

# Step 3: Union, then intersect with both datasets
candidate_genes = hvg_genes.union(marker_genes)
shared_genes = list(candidate_genes.intersection(sc_adata.var_names, visium_adata.var_names))
shared_genes = sorted(shared_genes)  # ensure deterministic order

# Step 4: Take top 2000
signature_genes = shared_genes[:n_genes_max]
print(f"Selected {len(signature_genes)} genes for Tangram.")

# Final subset
sc_adata = sc_adata[:, signature_genes].copy()
visium_adata = visium_adata[:, signature_genes].copy()

# -------------------------------
# Remove QC-filtered cells before dot plot
# -------------------------------
# Example: if QC-filtered cells are marked in the celltype_key column
qc_labels = ["QC_Filtered", "Doublets", "Low quality", 'Unknown III (SM)']  # adjust to your dataset
sc_adata = sc_adata[~sc_adata.obs[celltype_key].isin(qc_labels)].copy()

# -------------------------------
# Dot plot for signature gene expression per cell type
# -------------------------------

# # How many top genes to display per cell type in the dot plot
# n_top_plot = 5

# # Create dictionary: celltype -> top N marker genes (already computed above)
# top_marker_dict = {}
# for group in sc_adata.obs[celltype_key].unique():
#     df = sc.get.rank_genes_groups_df(sc_adata, group=group)
#     top_marker_dict[group] = df.head(n_top_plot).names.tolist()

# # Flatten gene list while keeping groupings
# gene_list = [g for genes in top_marker_dict.values() for g in genes if g in sc_adata.var_names]

# # Generate var_group_positions for grouping in the dot plot
# var_group_positions = []
# start_idx = 0
# for group, genes in top_marker_dict.items():
#     genes_present = [g for g in genes if g in sc_adata.var_names]
#     if genes_present:
#         end_idx = start_idx + len(genes_present) - 1
#         var_group_positions.append((start_idx, end_idx))
#         start_idx += len(genes_present)

# # Plot
# sc.pl.dotplot(
#     sc_adata,
#     var_names=gene_list,
#     groupby=celltype_key,
#     standard_scale="var",  # scale each gene to 0-1
#     dot_max=0.5,
#     dendrogram=False,
#     var_group_labels=list(top_marker_dict.keys()),
#     var_group_positions=var_group_positions
# )
# Build grouped gene dict using only genes present after filtering
n_top_plot = 5
grouped_var_names = {}
for group in sc_adata.obs[celltype_key].unique():
    df = sc.get.rank_genes_groups_df(sc_adata, group=group)
    genes_present = [g for g in df.head(n_top_plot).names.tolist() if g in sc_adata.var_names]
    if len(genes_present) > 0:
        grouped_var_names[group] = genes_present

# Assert we have at least one group and at least one gene
assert len(grouped_var_names) > 0, "No groups have marker genes after filtering"
assert sum(len(v) for v in grouped_var_names.values()) > 0, "No marker genes remain after filtering"

# Plot without manual group positions or labels
sc.pl.dotplot(
    sc_adata,
    var_names=grouped_var_names,      # dict triggers automatic grouping
    groupby=celltype_key,
    standard_scale="var",
    dot_max=0.5,
    dendrogram=False
)


# -------------------------------
# Run Tangram serially across batches
# -------------------------------
all_preds = []

for batch in visium_adata.obs[visium_batch_key].unique():
    print(f"\n[•] Processing batch: {batch}")
    output_path = os.path.join(output_dir, f"tangram_pred_{batch}.h5ad")
    
    if os.path.exists(output_path):
        print(f"[✓] Skipping batch {batch} (already computed)")
        batch_adata_map = sc.read_h5ad(output_path)
    else:
        batch_adata = visium_adata[visium_adata.obs[visium_batch_key] == batch].copy()
        sc_copy = sc_adata.copy()
        tg.pp_adatas(sc_copy, batch_adata, genes=signature_genes)
        batch_adata_map = tg.map_cells_to_space(sc_copy, batch_adata, mode='clusters', cluster_label=celltype_key, num_epochs = 500)
        batch_adata_map.write(output_path)
        print(f"[✔] Saved: {output_path}")

    all_preds.append(batch_adata_map)

# -------------------------------
# Merge final output
# -------------------------------
import anndata
merged_adata = anndata.concat(all_preds)
merged_adata.write_h5ad(os.path.join(output_dir, "tangram_prediction_merged.csv"))
print("\n✅ All Tangram batches completed and merged.")
os.makedirs('figures/tangram/', exist_ok = True)

import matplotlib.pyplot as plt
for batch in visium_adata.obs[visium_batch_key].unique():
    print(f"\n[•] Processing batch: {batch}")
    output_path = os.path.join(output_dir, f"tangram_pred_{batch}.h5ad")
    
    adata_st = visium_adata[visium_adata.obs[visium_batch_key] == batch].copy()

    if os.path.exists(output_path):
        ad_map = sc.read_h5ad(output_path)
    else:
        sc_copy = sc_adata.copy()
        tg.pp_adatas(sc_copy, batch_adata, genes=signature_genes)
        ad_map = tg.map_cells_to_space(sc_copy, batch_adata, mode='clusters', cluster_label=celltype_key, num_epochs = 500)
        ad_map.write(output_path)
        print(f"[✔] Saved: {output_path}")

    tg.project_cell_annotations(ad_map, adata_st, annotation=celltype_key)
    annotation_list = list(pd.unique(sc_adata.obs[celltype_key]))
    tg.plot_cell_annotation_sc(adata_st, annotation_list,perc=0.2, spot_size = 1, scale_factor = 1)

    plt.savefig(f'figures/tangram/{batch}.pdf', bbox_inches = 'tight')
    plt.savefig(f'figures/tangram/{batch}.png', bbox_inches = 'tight')
    plt.close()