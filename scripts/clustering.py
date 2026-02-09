import scanpy as sc
import scvi

adata = sc.read('results/scvi.h5ad')
adata = adata[:, ~(adata.var["mt"] | adata.var["ribo"]) | adata.var["hb"]]

sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=20, use_rep="X_scVI")
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.8)

adata.write('results/scvi.leiden.h5ad')