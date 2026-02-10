import scanpy as sc
import matplotlib.pyplot as plt

adata = sc.read('results/scanorama.leiden.grouped.h5ad')

tumor_map = {
    'PT01-1_NAT': 'Non-tumor',
    'PT01-2_TUM': 'Tumor',
    'PT01-3_NAT': 'Non-tumor',
    'PT01-4': 'Tumor',
    'PT01-5_NAT': 'Non-tumor',
    'PT01-6_TUM': 'Tumor',
    'PT05-2_TUM': 'Tumor',
    'PT05-4_TUM': 'Tumor',
}

adata.obs['tissue_type'] = adata.obs['batch'].str.replace('data/','').str.split('_binned_outputs').str[0].replace(tumor_map)

fig, ax = plt.subplots(1,1,figsize = (3,3))
sc.pl.umap(adata, color = ['meta_cluster'], save = 'umap_cluster.png', show = False, ax = ax)

fig, ax = plt.subplots(1,1,figsize = (3,3))
sc.pl.umap(adata, color = ['batch'], save = 'umap_batch.png', show = False, ax = ax)

fig, ax = plt.subplots(1,1,figsize = (3,3))
sc.pl.umap(adata, color = ['tissue_type'], save = 'umap_tissue_type.png', show = False, ax = ax)


sc.tl.rank_genes_groups(adata, 'meta_cluster', method='wilcoxon', use_raw = False)
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, show = False)
plt.savefig('figures/meta_cluster_signature_genes.pdf')
plt.tight_layout()
plt.close()

marker_df = sc.get.rank_genes_groups_df(adata, group=None)
print(marker_df.groupby('group').head(3))

idx = adata.obs['meta_cluster'].value_counts() > 1000
adata_subset = adata[adata.obs['meta_cluster'].isin(idx[idx].index)].copy()
adata_subset.obs['meta_cluster'] = pd.Categorical(adata_subset.obs['meta_cluster'], categories = idx[idx].index.tolist())


sc.pl.rank_genes_groups_heatmap(
    adata,
    groupby='meta_cluster',
    n_genes=10,
    use_raw=False,
    show_gene_labels=True,
    swap_axes=True,
    show = False
)
plt.savefig('figures/genes_meta_cluster_heatmap.pdf')
plt.tight_layout()
plt.close()

# Extract top 5 marker genes per cluster
markers_df = sc.get.rank_genes_groups_df(adata, group=None)

top_genes = (
    markers_df.groupby("group")
    .apply(lambda x: x.nlargest(5, "logfoldchanges"))
    .names.unique()
    .tolist()
)
marker_dict = top_genes.groupby('group')['names'].apply(list).to_dict()

# Plot heatmap of selected genes
sc.pl.heatmap(
    adata,
    var_names=marker_dict,
    groupby="meta_cluster",
    use_raw=False,
    swap_axes=True,
    show_gene_labels=True,
    standard_scale='var',  # z-score per gene
    show = False
)
plt.savefig('figures/top_genes_meta_cluster_heatmap.pdf')
plt.tight_layout()
plt.close()

sc.pl.tracksplot(
    adata,
    var_names=marker_dict,
    groupby="meta_cluster",
    show = False,
    dendrogram=False,
    sort=True
)
plt.savefig('figures/top_genes_meta_cluster_tracksplot.pdf')
plt.tight_layout()
plt.close()

sc.pl.dotplot(adata, var_names=marker_dict, groupby="meta_cluster", standard_scale='var', show = False)
plt.savefig('figures/top_genes_meta_cluster_dotplot.pdf')
plt.tight_layout()
plt.close()

sc.pl.matrixplot(adata, var_names=marker_dict, groupby="meta_cluster", standard_scale='var', show = False)
plt.savefig('figures/top_genes_meta_cluster_matrixplot.pdf')
plt.tight_layout()
plt.close()

sc.pl.spatial(adata, color=['meta_cluster'], spot_size=1.5)

# Compute spatial neighbors
sq.gr.spatial_neighbors(adata)
sq.gr.spatial_autocorr(adata, mode='moran')
sq.pl.spatial_autocorr(adata, mode='moran', genes=adata.var_names[:4])
