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
sc_data_file = 'output/sc_reference.h5ad'
visium_hd_data_file = 'output/adata_filtered.h5ad'
celltype_key = 'Level1'
visium_batch_key = 'sample'
sc_batch_key = 'BC'
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
    tg.plot_cell_annotation_sc(adata_st, annotation_list,perc=0.2, spot_size = 20, scale_factor = 1)

    plt.savefig(f'figures/tangram/{batch}.pdf', bbox_inches = 'tight')
    plt.savefig(f'figures/tangram/{batch}.png', bbox_inches = 'tight')
    plt.close()