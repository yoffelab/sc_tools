import scanpy as sc
from sklearn.cluster import KMeans
import pandas as pd

def group_small_clusters(
    adata,
    cluster_key="leiden",
    min_cells=1000,
    small_label="Other",
    new_key=None,
    inplace=True,
):
    """
    Collapse rare clusters into `small_label` and order categories by size.

    Parameters
    ----------
    adata : AnnData
        Object containing clustering results.
    cluster_key : str
        Column in `adata.obs` holding cluster labels.
    min_cells : int
        Threshold below which clusters are considered small.
    small_label : str
        Label assigned to merged small clusters.
    new_key : str or None
        Name of the output column. If None, appends '_grouped' to `cluster_key`.
    inplace : bool
        If True, write the new column to `adata.obs`; otherwise return a Series.

    Returns
    -------
    None or pandas.Series
        When `inplace=True`, the function edits `adata` in place and returns None.
        Otherwise it returns the new cluster assignments as a Series.
    """
    if new_key is None:
        new_key = f"{cluster_key}_grouped"

    # Current cluster labels as strings
    original = adata.obs[cluster_key].astype(str)

    # Identify small clusters
    counts = original.value_counts()
    small_clusters = counts[counts < min_cells].index

    # Vectorised reassignment
    grouped = original.where(~original.isin(small_clusters), other=small_label)

    # Order categories by size (descending)
    ordered_cats = (
        grouped.value_counts(sort=True)
        .sort_values(ascending=False)
        .index
        .tolist()
    )
    grouped = pd.Categorical(grouped, categories=ordered_cats, ordered=True)

    if inplace:
        adata.obs[new_key] = grouped
    else:
        return grouped

# import numpy as np
# import pandas as pd
# from sklearn.metrics import pairwise_distances
# from matplotlib.colors import ListedColormap, to_hex
# import matplotlib.pyplot as plt
# import seaborn as sns
# import re
import re
import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances
from sklearn.cluster import AgglomerativeClustering
import colorspacious as cs
import matplotlib.colors as mcolors
import seaborn as sns

def assign_cluster_colors(
    adata,
    cluster_key="leiden",
    embed_key="X_pca",
    n_dims=30,
    n_hues=8,
    metric="cosine",
    color_key=None,
    seed=42,
    min_per_class=3,
    inplace=True,
):
    """
    Assign interpretable, colour-blind-safe colours to clusters or cell types.

    Strategy:
    1. If labels match known cell types (e.g. CD8, fibro), use interpretable colours.
    2. Otherwise, use fallback: distance-aware colour-blind-safe palette.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    cluster_key : str
        Column in adata.obs containing cluster or cell type labels.
    embed_key : str
        Key in adata.obsm for the reduced coordinates (e.g., X_pca).
    n_dims : int
        Number of dimensions to use from embed_key.
    n_hues : int
        Number of colour-blind-safe base hues (max 8 for Okabe–Ito palette).
    metric : str
        Distance metric for transcriptomic distance (e.g., "cosine", "euclidean").
    color_key : str or None
        Where to store colours in adata.uns (default: {cluster_key}_colors).
    seed : int
        Random seed for reproducibility.
    min_per_class : int
        Minimum shades per biological class in interpretable mode.
    inplace : bool
        Whether to store the result in adata.uns.

    Returns
    -------
    None or dict
        Palette stored in adata.uns or returned as dictionary if inplace=False.
    """
    if color_key is None:
        color_key = f"{cluster_key}_colors"

    rng = np.random.default_rng(seed)
    labels = adata.obs[cluster_key].astype(str)
    categories = labels.unique()

    # ──────────────────────────────────────────────────────────────────────
    # Step 1: Determine if any known biological keyword is present
    # ──────────────────────────────────────────────────────────────────────
    class_patterns = {
        "epithelial": r"(epi|tumou?r|alveo|club|acinar)",
        "stromal": r"(strom|fibro|mesench|caf|endothel)",
        "lymphoid": r"(cd4|cd8|treg|b[\s_\-]?cell|lymph|t[_\s\-]?|nk)",
        "myeloid": r"(macro|mono|dend|mast|neutro|myeloid|microglia)",
    }

    match_found = any(re.search("|".join(class_patterns.values()), c.lower()) for c in categories)

    # ──────────────────────────────────────────────────────────────────────
    # Case A: Use interpretable palette if cell-type keywords are present
    # ──────────────────────────────────────────────────────────────────────
    if match_found:
        coords_df = pd.DataFrame(adata.obsm[embed_key][:, :n_dims], index=labels)
        centroids = coords_df.groupby(level=0).mean()
        
        # Classify each cluster into a category
        class_map = {}
        for cat in categories:
            for cls, pattern in class_patterns.items():
                if re.search(pattern, cat.lower()):
                    class_map[cat] = cls
                    break
            else:
                class_map[cat] = "other"

        # Predefined color palettes per class
        class_palettes = {
            "epithelial": sns.color_palette("Reds", min_per_class),
            "stromal": sns.color_palette("Greens", min_per_class),
            "lymphoid": sns.color_palette("Blues", min_per_class),
            "myeloid": sns.color_palette("Purples", min_per_class),
            "other": sns.color_palette("gray", min_per_class),
        }

        palette = {}

        # Assign colours within each class based on distance to centroid
        for cls, hue_list in class_palettes.items():
            members = [c for c, assigned in class_map.items() if assigned == cls]
            if not members:
                continue

            if len(members) == 1:
                palette[members[0]] = mcolors.to_hex(hue_list[0])
                continue

            sub_coords = centroids.loc[members].to_numpy()
            dist = pairwise_distances(sub_coords, sub_coords.mean(0, keepdims=True), metric=metric).ravel()
            sorted_indices = np.argsort(dist)
            sorted_members = [members[i] for i in sorted_indices]
            sorted_colors = sns.color_palette(hue_list, len(sorted_members))

            for cat, col in zip(sorted_members, sorted_colors):
                palette[cat] = mcolors.to_hex(col)

    # ──────────────────────────────────────────────────────────────────────
    # Case B: Use fallback distance-aware CB-safe palette
    # ──────────────────────────────────────────────────────────────────────
    else:
        coords = adata.obsm[embed_key][:, :n_dims]
        centroid_df = pd.DataFrame(coords, index=labels)
        centroids = centroid_df.groupby(level=0).mean().to_numpy()

        if len(categories) <= n_hues:
            hue_ids = np.arange(len(categories))
        else:
            cluster = AgglomerativeClustering(n_clusters=n_hues, affinity=metric, linkage="average")
            hue_ids = cluster.fit_predict(centroids)

        okabe_ito_base = np.array([
            "#E69F00", "#56B4E9", "#009E73", "#F0E442",
            "#0072B2", "#D55E00", "#CC79A7", "#999999"
        ])[:n_hues]

        palette = {}

        for h in range(n_hues):
            members = np.where(hue_ids == h)[0]
            if not len(members):
                continue

            sub_centroids = centroids[members]
            centre = sub_centroids.mean(0, keepdims=True)
            dist = pairwise_distances(sub_centroids, centre, metric=metric).ravel()
            sorted_indices = members[np.argsort(dist)]

            # Convert base hue to OKLab
            base_rgb = np.array(mcolors.to_rgb(okabe_ito_base[h]))[None, :]
            base_oklab = cs.cspace_convert(base_rgb, "sRGB1", "OKLab").ravel()
            L_range = np.linspace(0.95, 0.25, len(sorted_indices))
            a, b = base_oklab[1], base_oklab[2]
            oklab = np.column_stack([L_range, np.full(len(L_range), a), np.full(len(L_range), b)])
            rgb = cs.cspace_convert(oklab, "OKLab", "sRGB1")
            rgb = np.clip(rgb, 0, 1)
            hex_colors = [mcolors.to_hex(c) for c in rgb]

            for idx, hex_col in zip(sorted_indices, hex_colors):
                palette[categories[idx]] = hex_col

    # ──────────────────────────────────────────────────────────────────────
    # Final output
    # ──────────────────────────────────────────────────────────────────────
    sorted_palette = [palette[c] for c in sorted(categories)]
    
    if inplace:
        adata.uns[color_key] = sorted_palette
    else:
        return {c: palette[c] for c in sorted(categories)}



adata = sc.read('output/adata_corrected_scvi.h5ad')
for res in [0.001, 0.005, 0.01]:
    sc.tl.leiden(adata, resolution = res, key_added = f'leiden_{res}')
    
for res in [0.001, 0.005, 0.01]:
    group_small_clusters(adata, cluster_key = f'leiden_{res}', min_cells = 10000)
    assign_cluster_colors(
        adata,
        cluster_key=f'leiden_{res}_grouped',
        embed_key="X_scVI",
        n_dims=5,
        min_per_class=4,
    )
    sc.pl.umap(adata, color = [f'leiden_{res}_grouped'], show = False, save = f'umap_leiden_{res}.png', legend_loc="on data",
    legend_fontsize=10,
    legend_fontoutline=2,)

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

adata.obs['tissue_type'] = adata.obs['sample'].str.split('_binned_outputs').str[0].replace(tumor_map)
sc.pl.umap(adata, color = ['sample'], show = False, save = 'umap_sample.png')

adata.write('output/scvi.leiden.h5ad')

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
import seaborn as sns
import matplotlib.pyplot as plt
adata = sc.read('output/scvi.leiden.h5ad')
cluster_labels = adata.obs['leiden']
means = adata.to_df().groupby(cluster_labels).mean()

# Hierarchical clustering on cluster means
Z = linkage(means, method='ward')
# dendrogram(Z, labels=means.index)
# plt.show()

# Cut dendrogram at k groups
k = 15
meta_labels = fcluster(Z, k, criterion='maxclust')
cluster_to_meta = dict(zip(means.index, meta_labels.astype(str)))

adata.obs['meta_cluster'] = adata.obs['leiden'].map(cluster_to_meta)
adata.obs['meta_cluster'] = adata.obs['meta_cluster'].astype('category')
adata.write('output/scvi.leiden.grouped.h5ad')

