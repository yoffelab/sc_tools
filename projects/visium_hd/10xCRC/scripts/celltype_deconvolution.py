import scanpy as sc
import scvi
import cell2location
import tangram as tg
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

# -------------------------------
# Configuration
# -------------------------------
sc_data_file = 'output/sc_reference.h5ad'
visium_hd_data_file = 'output/adata_filtered.h5ad'

celltype_key = 'Level1'  # Could also be 'Level2'
batch_key = 'BC'

# -------------------------------
# Step 1: Load and Prepare scRNA-seq Data
# -------------------------------
sc_adata = sc.read_h5ad(sc_data_file)
sc_adata.var_names_make_unique()

# Set up for scVI (batch-aware, label-aware)
scvi_dir = 'models/scvi_model/'
if not os.path.exists(scvi_dir + 'model.pt'):
    scvi.model.SCVI.setup_anndata(sc_adata, batch_key=batch_key, labels_key=celltype_key)
    sc_model = scvi.model.SCVI(sc_adata)
    sc_model.train()
    sc_model.save(scvi_dir, overwrite=True)
else:
    sc_model = scvi.model.SCVI.load(scvi_dir, adata=sc_adata)

# -------------------------------
# Step 2: Train Cell2location Reference Model
# -------------------------------
from cell2location.models import RegressionModel

ref_model_dir = 'models/cell2location_reference/'

sc_adata_copy = sc_adata.copy()
RegressionModel.setup_anndata(sc_adata_copy, labels_key=celltype_key)

if not os.path.exists(ref_model_dir + 'model.pt'):
    ref_model = RegressionModel(sc_adata_copy)
    ref_model.train(max_epochs=250)
    ref_model.save(ref_model_dir, overwrite=True)
    sc_adata_copy = ref_model.export_posterior(sc_adata_copy)
    sc_adata_copy.write('output/sc_adata_reference_post.h5ad')
else:
    ref_model = RegressionModel.load(ref_model_dir, adata=sc_adata_copy)
    sc_adata_copy = sc.read_h5ad('output/sc_adata_reference_post.h5ad')

# -------------------------------
# Step 3: Load Visium HD Data
# -------------------------------
visium_adata = sc.read_h5ad(visium_hd_data_file)
visium_adata.var_names_make_unique()
visium_adata.obs_names_make_unique()

# If raw exists, use its counts but match current gene space
if visium_adata.raw is not None:
    raw_adata = visium_adata.raw.to_adata()
    visium_adata.X = raw_adata[:, visium_adata.var_names].X

# Subset to shared genes
shared_genes = sc_adata.var_names.intersection(visium_adata.var_names)
visium_adata = visium_adata[:, shared_genes]

# Sanity check
if hasattr(visium_adata.X, "data"):
    assert np.allclose(visium_adata.X.data.astype(int), visium_adata.X.data), "adata.X contains non-integers"
else:
    assert np.allclose(visium_adata.X.astype(int), visium_adata.X), "adata.X contains non-integers"

# -------------------------------
# Step 4: Run Cell2location Mapping
# -------------------------------

c2l_model_dir = 'models/cell2location_mapping/'

cell2location.models.Cell2location.setup_anndata(visium_adata)

cell_state_df = sc_adata_copy.varm["means_per_cluster_mu_fg"]
cell_state_df = pd.DataFrame(cell_state_df, index=sc_adata_copy.var_names)
cell_state_df = cell_state_df.loc[visium_adata.var_names]

if not os.path.exists(c2l_model_dir + 'model.pt'):
    c2l_model = cell2location.models.Cell2location(
        visium_adata,
        cell_state_df=cell_state_df
    )
    c2l_model.train(max_epochs=5000, batch_size=128, train_size=1)
    c2l_model.save(c2l_model_dir, overwrite=True)
    visium_adata = c2l_model.export_posterior(visium_adata)
    visium_adata.write('output/visium_cell2location_post.h5ad')
else:
    c2l_model = cell2location.models.Cell2location.load(c2l_model_dir, adata=visium_adata)
    visium_adata = sc.read_h5ad('output/visium_cell2location_post.h5ad')

# -------------------------------
# Step 5: Visualize Cell Type Abundances
# -------------------------------
cell_types = visium_adata.obsm['q05_cell_abundance_w_sf'].columns

for cell in cell_types:
    sc.pl.spatial(
        visium_adata,
        color=[cell],
        spot_size=1.2,
        cmap="magma",
        title=f"{cell} abundance"
    )

# -------------------------------
# Step 6: Tangram (Alternative Mapping)
# -------------------------------
# Ensure shared genes are registered to both AnnData objects
tg.pp_adatas(sc_adata, visium_adata, genes=shared_genes)

# Run Tangram in 'clusters' mode, using the provided cell type label
tg.map_cells_to_space(sc_adata, visium_adata, mode='clusters', cluster_label=celltype_key)
