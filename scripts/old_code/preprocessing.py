import scanpy as sc
from sklearn.cluster import KMeans
import pandas as pd

# adata = sc.read('results/scanorama_corrected.h5ad')
adata = sc.read('results/harmony.h5ad')

# sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
sc.pl.umap(adata, color = ['leiden','batch'], show = False, save = 'umap')

adata.write('results/scanorama.leiden.h5ad')


from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
import seaborn as sns
import matplotlib.pyplot as plt
adata = sc.read('results/scanorama.leiden.h5ad')
cluster_labels = adata.obs['leiden']
means = pd.DataFrame(
    adata.X,
    index=adata.obs_names,
    columns=adata.var_names
).groupby(cluster_labels).mean()

# Hierarchical clustering on cluster means
Z = linkage(means, method='ward')
dendrogram(Z, labels=means.index)
plt.show()

# Cut dendrogram at k groups
k = 15
meta_labels = fcluster(Z, k, criterion='maxclust')
cluster_to_meta = dict(zip(means.index, meta_labels.astype(str)))

adata.obs['meta_cluster'] = adata.obs['leiden'].map(cluster_to_meta)
adata.obs['meta_cluster'] = adata.obs['meta_cluster'].astype('category')
adata.write('results/scanorama.leiden.grouped.h5ad')
