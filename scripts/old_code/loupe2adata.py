import scanpy as sc
import squidpy as sq
import pandas as pd
import anndata
from glob import glob
from tqdm import tqdm
import numpy as np
import os
import sthd

glob_pattern = 'data/*/outs/binned_outputs/square_008um/'
loupe_files = sorted(glob(glob_pattern))
location_files = sorted(glob(f'{glob_pattern}/spatial/tissue_positions.parquet'))

raw_filename = 'results/raw.008um.concat.h5ad'
if not os.path.exists(raw_filename):
    adatas = list()
    for path, loc_path in tqdm(zip(loupe_files,location_files)):

        a_ = sq.read.visium(path, load_images = False, counts_file = 'filtered_feature_bc_matrix.h5')
        
        df = pd.read_parquet(loc_path)
        df = df.set_index('barcode')
        a_.obs = df.loc[a_.obs.index]

        a_.obs['batch'] = path
        a_.var_names_make_unique()

        adatas.append(a_)

    adata = anndata.concat(adatas, merge = 'same', uns_merge = 'same')
    adata.obsm['spatial'] = np.array(adata.obs[['pxl_col_in_fullres', 'pxl_row_in_fullres']])

    adata.raw = adata.copy()
    adata.write(raw_filename)




adata_raw = sc.read(raw_filename)
sthd.qc.calculate_qc_metrics(adata_raw)
sthd.pl.qc_overview(adata_raw, save_name="qc_pre_filter.pdf")
sthd.io.save_h5ad(adata_raw, "output/raw_data.h5ad")

adata_raw = sc.read("output/raw_data.h5ad")

# process metadata
meta = pd.read_excel("metadata/sample_metadata.xlsx")
obs = adata_raw.obs.copy()
obs["sample"] = (
    obs["batch"].str.split("/").str[1].str.replace(r"_S\d$", "", regex=True)
)

# 3) Join without changing obs index
meta_idx = meta.set_index("ROBIN Sample info")
obs = obs.join(meta_idx, on="sample", how="left", rsuffix="_meta")
adata_raw.obs = obs

adata_raw.obs['sample_id'] = adata_raw.obs['Batch'] + \
    '_P' + adata_raw.obs['Patient#'].astype(str) + \
    '_' + adata_raw.obs['Timepoint'] + \
    '_' + adata_raw.obs['Tumor'].astype(str) + \
    '_' + adata_raw.obs['Organ'] + \
    '_' + adata_raw.obs['Fixation']



high_res_image_files = sorted(glob('data/*/outs/spatial/tissue_hires_image.png'))
low_res_images_files = sorted(glob('data/*/outs/spatial/tissue_lowres_image.png'))
'''
from PIL import Image
adata_raw.obs['library_id'] = adata_raw.obs['batch'].str.split('/').str[1]
spatial_dict = dict()
for library_id in tqdm(adata_raw.obs['library_id'].unique()):
    
    spatial_dict[library_id] ={
        "images": {
            "hires": Image.open(f'data/{library_id}/outs/spatial/tissue_hires_image.png'),
            "lowres": Image.open(f'data/{library_id}/outs/spatial/tissue_lowres_image.png')
        },
        "scalefactors": 1,
        "metadata": {}
    }

adata_raw.uns['spatial'] = spatial_dict

if 'spatial' in adata.uns:
    del adata.uns['spatial']
adata.obs['library_id'] = adata.obs['batch'].str.split('/').str[1]
spatial_dict = dict()
for library_id in tqdm(adata.obs['library_id'].unique()):
    
    sf_p = f'data/{library_id}/outs/binned_outputs/square_008um/spatial/scalefactors_json.json'
    with open(sf_p, "r") as f:
        sf = json.load(f)
        # required keys for Scanpy
    assert "tissue_hires_scalef" in sf and "tissue_lowres_scalef" in sf, \
        f"Missing hires or lowres scalefactors in {sf_p}"

    spatial_dict[library_id] = {
        "images": {
            "hires": np.array(Image.open(f'data/{library_id}/outs/spatial/tissue_hires_image.png')),
            "lowres": np.array(Image.open(f'data/{library_id}/outs/spatial/tissue_lowres_image.png')
        },
        "scalefactors": {
            "tissue_hires_scalef": float(sf["tissue_hires_scalef"]),
            "tissue_lowres_scalef": float(sf["tissue_lowres_scalef"]),
            # keep if present, else reasonable fallback
            "spot_diameter_fullres": float(sf.get("spot_diameter_fullres", 50.0)),
        },
        "metadata": {}
    }

adata.uns['spatial'] = spatial_dict


for library_id in tqdm(adata.obs['library_id'].unique()[:2]):
    a_ = adata[adata.obs['library_id'] == library_id].copy()
    sc.pl.spatial(a_, color = 'leiden', library_id = library_id)

sc.pl.spatial(adata, color = 'leiden', library_id = 'PT01-2_TUM')

# TODO
# 1. test whether spatial plotting works
# 2. save images to all AnnData
'''

# QC plot
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
# Replace 'sample_id' with your actual batch key if different
df = adata_raw.obs[['log1p_total_counts', 'sample_id']].copy()
df = df.rename(columns={'log1p_total_counts': 'log_counts'})

# Optional: Sort sample_ids
df['sample_id'] = pd.Categorical(df['sample_id'], ordered=True)

# Create ridgeplot (joyplot)
g = sns.FacetGrid(
    df,
    row='sample_id',
    hue='sample_id',
    palette = 'Set2',
    aspect=4,
    height=1.2,
    sharey=False,
)

g.map(sns.kdeplot, 'log_counts', fill=True, alpha=0.7, bw_adjust=1)
g.set_titles(row_template='{row_name}')
g.set(xlim=(0, None), xlabel='log1p(total counts)', ylabel='')
g.despine(left=True)

# Natural values and their log1p
raw_vals = [10, 20, 50, 100]
log_vals = [np.log1p(v) for v in raw_vals]

# Loop over all axes in the FacetGrid
for ax in g.axes.flatten():
    for lv, rv in zip(log_vals, raw_vals):
        ax.axvline(x=lv, color='gray', linestyle='dotted', linewidth=1)
        ax.text(lv, ax.get_ylim()[1]*0.5, f'{rv}', rotation=90, ha='right', va='top', fontsize=7, color='gray')

plt.tight_layout()
plt.savefig('figures/QC per sample.pdf', bbox_inches = 'tight')
plt.close()



# Filter samples that failed QC
adata_raw = adata_raw[adata_raw.obs['QC'] != 'Failed'].copy()

# Filtering & normalization
adata_filtered = sthd.pp.filter_cells_and_genes(adata_raw)
sthd.pp.normalize_and_log(adata_filtered)
sthd.qc.calculate_qc_metrics(adata_filtered)
sthd.pl.qc_overview(adata_filtered, save_name="qc_post_filter.pdf")

if 'QC' in adata_filtered.obs:
    del adata_filtered.obs['QC']
sthd.io.save_h5ad(adata_filtered, "output/adata_filtered.h5ad")

# HVG & SVG
adata_filtered = sthd.io.read_h5ad("output/adata_filtered.h5ad")
adata_filtered = sthd.pp.filter_by_consensus_hvg_svg(adata_filtered, hvg_top_n=5000, svg_top_n=5000,)
# filter mt, rb, and, hb genes
adata = adata[:, ~(adata.var["mt"] | adata.var["ribo"]) | adata.var["hb"]]
sthd.io.save_h5ad(adata_filtered, "output/adata_vg_filtered.h5ad")


adata_batch = sthd.pp.batch_integrate(adata_filtered, batch_key = 'batch', method = 'scvi')
adata_batch.write('results/scvi.h5ad')
