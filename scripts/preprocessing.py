
import numpy as np
import pandas as pd
import scanpy as sc
import sthd

# adata_raw = sc.read(raw_filename)
# sthd.qc.calculate_qc_metrics(adata_raw)
# sthd.pl.qc_overview(adata_raw, save_name="qc_pre_filter.pdf")
# sthd.io.save_h5ad(adata_raw, "output/raw_data.h5ad")

adata = sc.read("results/adata.annotated.p2.h5ad")

# process metadata

# 3) Join without changing obs index
# meta_idx = meta.set_index("ROBIN Sample info")
# obs = obs.join(meta_idx, on="sample", how="left", rsuffix="_meta")
# adata_raw.obs = obs

# TODO
# 1. test whether spatial plotting works
# 2. save images to all AnnData

# QC plot
import matplotlib.pyplot as plt
import seaborn as sns

adata_raw = adata
# Replace 'sample_id' with your actual batch key if different
df = adata_raw.obs[["log1p_total_counts", "library_id"]].copy()
df = df.rename(columns={"log1p_total_counts": "log_counts"})

# Optional: Sort library_ids
df["library_id"] = pd.Categorical(df["library_id"], ordered=True)

# Create ridgeplot (joyplot)
g = sns.FacetGrid(
    df,
    row="library_id",
    hue="library_id",
    palette="Set2",
    aspect=4,
    height=1.2,
    sharey=False,
)

g.map(sns.kdeplot, "log_counts", fill=True, alpha=0.7, bw_adjust=1)
g.set_titles(row_template="{row_name}")
g.set(xlim=(0, None), xlabel="log1p(total counts)", ylabel="")
g.despine(left=True)

# Natural values and their log1p
raw_vals = [10, 20, 50, 100]
log_vals = [np.log1p(v) for v in raw_vals]

# Loop over all axes in the FacetGrid
for ax in g.axes.flatten():
    for lv, rv in zip(log_vals, raw_vals):
        ax.axvline(x=lv, color="gray", linestyle="dotted", linewidth=1)
        ax.text(
            lv,
            ax.get_ylim()[1] * 0.5,
            f"{rv}",
            rotation=90,
            ha="right",
            va="top",
            fontsize=7,
            color="gray",
        )

plt.tight_layout()
plt.savefig("figures/QC per sample.pdf", bbox_inches="tight")
plt.close()


# Filtering & normalization
adata_filtered = sthd.pp.filter_cells_and_genes(adata_raw)
sthd.pp.normalize_and_log(adata_filtered)
sthd.qc.calculate_qc_metrics(adata_filtered)
sthd.pl.qc_overview(adata_filtered, save_name="qc_post_filter.pdf")

if "QC" in adata_filtered.obs:
    del adata_filtered.obs["QC"]

adata_filtered.obs["sample"] = adata_filtered.obs["library_id"]
sthd.io.save_h5ad(adata_filtered, "results/adata_filtered.h5ad")

# HVG & SVG
adata_filtered = sthd.io.read_h5ad("results/adata_filtered.h5ad")
adata_filtered = sthd.pp.filter_by_consensus_hvg_svg(
    adata_filtered,
    hvg_top_n=5000,
    svg_top_n=5000,
)
# filter mt, rb, and, hb genes
adata = adata[:, ~(adata.var["mt"] | adata.var["ribo"]) | adata.var["hb"]]
sthd.io.save_h5ad(adata_filtered, "results/adata.normalized.p3.h5ad")


adata_batch = sthd.pp.batch_integrate(adata_filtered, batch_key="library_id", method="scvi")
adata_batch.write("results/scvi.h5ad")
