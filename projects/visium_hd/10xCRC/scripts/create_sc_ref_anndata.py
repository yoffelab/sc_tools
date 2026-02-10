import scanpy as sc
import pandas as pd
import anndata
import numpy as np

adata = sc.read_10x_h5("data/single_cell/HumanColonCancer_VisiumHD-main/AggrOutput/outs/count/filtered_feature_bc_matrix.h5")
metadata = pd.read_csv('data/single_cell/HumanColonCancer_VisiumHD-main/Metadata/SingleCell_MetaData.csv', index_col = 0)

# Optionally: make var names unique and inspect
adata.var_names_make_unique()

adata.obs = metadata
adata.obsm['X_umap'] = np.array(adata.obs[['UMAP1', 'UMAP2']])
adata.write('output/sc_reference.h5ad')