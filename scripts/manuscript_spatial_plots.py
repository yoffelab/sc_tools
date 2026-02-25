import json
import os

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from tqdm import tqdm

sc.settings.set_figure_params(dpi=300, dpi_save=400)

# Assuming 'data.json' contains valid JSON data
with open("metadata/gene_signatures.json") as file:
    spatial_signatures = json.load(file)

adata = sc.read("results/adata.normalized.scored.p35.h5ad")


# Spatial Plots
os.makedirs("figures/manuscript/EMT_spatial/", exist_ok=True)
for lib in tqdm(adata.obs["library_id"].unique()):
    a = adata[adata.obs["library_id"] == lib].copy()
    sc.pl.spatial(
        a,
        color=spatial_signatures["Tumor_Cells"]["EMT_Tumor"],
        frameon=False,
        library_id=lib,
        cmap="coolwarm",
        show=False,
    )
    plt.savefig(f"figures/manuscript/EMT_spatial/{lib}.pdf", bbox_inches="tight")
    plt.savefig(f"figures/manuscript/EMT_spatial/{lib}.png", bbox_inches="tight")
    plt.close()


# EMT Related Genes
emt_dict = {
    "AT2 Tumor": spatial_signatures["Tumor_Cells"]["Alveolar_Type2_Tumor"],
    "AT1 Tumor": spatial_signatures["Tumor_Cells"]["Alveolar_Type1_Tumor"],
    "Proliferative Tumor": spatial_signatures["Tumor_Cells"]["Proliferative_Tumor"],
    "EMT": spatial_signatures["Tumor_Cells"]["EMT_Tumor"],
}

for k in emt_dict.keys():
    emt_dict[k] = [gene for gene in emt_dict[k] if gene in adata.var.index]

pa = {
    "Solid": "Solid Tumor",
    "Non-Solid": "Non-Solid Tumor",
    "Normal": "Normal Alveolar Cells",
    "Solid Blood Vessel": "Solid Blood Vessel",
    "Solid Bronchus": "Solid Bronchus",
    "Solid Scar Tissue": "Solid Scar Tissue",
    "Non-Solid Blood Vessel": "Non-Solid Blood Vessel",
    "Non-Solid Bronchus": "Non-Solid Bronchus",
    "Normal Blood Vessel": "Normal Blood Vessel",
    "Normal Bronchus": "Normal Bronchus",
    "Solid TLS": "TLS Solid",
    "TLS Solid": "TLS Solid",
    "Non-Solid TLS": "TLS Non-Solid",
    "TLS Non-Solid": "TLS Non-Solid",
    "TLS Normal": "TLS Normal",
}


def keep_first_unique(input_list):
    seen = set()
    unique_list = []
    for item in input_list:
        if item not in seen:
            unique_list.append(item)
            seen.add(item)
    return unique_list


adata.obs["pa"] = pd.Categorical(
    adata.obs["pathologist_annotation"].astype(str).replace(pa),
    categories=keep_first_unique(pa.values()),
)

sc.pl.matrixplot(
    adata,
    groupby="pa",  # ['solidity_type', 'architecture_type'],
    var_names=emt_dict,
    var_group_rotation=0,
    cmap="coolwarm",
    standard_scale="var",
    show=False,
)
plt.savefig("figures/manuscript/EMT_architecture_markers_matrixplot.pdf", bbox_inches="tight")
plt.savefig("figures/manuscript/EMT_architecture_markers_matrixplot.png", bbox_inches="tight")
plt.close()

adata

import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns
from tqdm import tqdm

adata.obs["Tumor Score"] = adata.obs[
    [
        "sig:Tumor_Cells/Alveolar_Type2_Tumor_z",
        "sig:Tumor_Cells/Alveolar_Type1_Tumor_z",
        "sig:Tumor_Cells/Proliferative_Tumor_z",
    ]
].mean(axis=1)
adata.obs["EMT Score"] = adata.obs["sig:Tumor_Cells/EMT_Tumor_z"]

categories = ["Normal Alveolar Cells", "Non-Solid Tumor", "Solid Tumor"]
a = adata[adata.obs["pa"].isin(categories)].copy()
a.obs["pa"] = pd.Categorical(a.obs["pa"], categories=categories)

# EMT Related Gene Signatures
keys = [
    "sig:Tumor_Cells/Alveolar_Type2_Tumor_z",
    "sig:Tumor_Cells/Alveolar_Type1_Tumor_z",
    "sig:Tumor_Cells/Proliferative_Tumor_z",
    "sig:Tumor_Cells/EMT_Tumor_z",
]

ylabels = ["AT2 Tumor score", "AT1 Tumor score", "Prol. Tumor score", "EMT Tumor score"]
i = 0
fig, axes = plt.subplots(nrows=1, ncols=4, figsize=(9, 2.5))
for ax, key in zip(axes, keys):
    sc.pl.violin(a, keys=key, groupby="pa", show=False, ax=ax, legend=None, jitter=0.1, size=0.3)
    ax.set_ylabel(ylabels[i])
    ax.set_xlabel("")
    ax.set_xticklabels(["N", "NS", "S"])
    for coll in ax.collections:
        if hasattr(coll, "get_offsets"):  # only scatter objects
            coll.set_alpha(0.6)

    i += 1
plt.tight_layout()
plt.savefig("figures/manuscript/EMT tumor scores.png", bbox_inches="tight")
plt.savefig("figures/manuscript/EMT tumor scores.pdf", bbox_inches="tight")
plt.close()


for lib in tqdm(a.obs["library_id"].unique()):
    a_ = a[a.obs["library_id"] == lib].copy()
    sc.pl.spatial(a_, color=keys, frameon=False, library_id=lib, cmap="coolwarm", show=False)
    plt.savefig(f"figures/manuscript/EMT_spatial/{lib}_scores.pdf", bbox_inches="tight")
    plt.savefig(f"figures/manuscript/EMT_spatial/{lib}_scores.png", bbox_inches="tight")
    plt.close()


fig, axes = plt.subplots(1, 3, dpi=300, figsize=(6, 2.5))
for i, ax in tqdm(enumerate(axes.flatten())):
    a_tmp = a[a.obs["pa"] == categories[i]].copy()
    sns.kdeplot(a_tmp.obs, x="Tumor Score", y="EMT Score", fill=True, ax=ax, cmap="YlOrRd")
    ax.set_title(categories[i])
    ax.set_xlim(-2.6, 5.2)
    ax.set_ylim(-2.6, 5.2)
    ax.grid(False)
    if i > 0:
        ax.set_ylabel("")
    # ax.legend().remove()
    sns.despine()
plt.tight_layout()
plt.savefig("figures/manuscript/EMT_architecture_score_kdeplot.png", bbox_inches="tight")
plt.savefig("figures/manuscript/EMT_architecture_score_kdeplot.pdf", bbox_inches="tight")
plt.close()
