# in conda env robin / sthd

import scanpy as sc
import squidpy as sq
import pandas as pd
import anndata
from glob import glob
from tqdm import tqdm
import numpy as np
import os
import sthd

import json
from PIL import Image
import matplotlib.pyplot as plt
from pathlib import Path
import re
sc.settings.set_figure_params(dpi=300, dpi_save = 300) 

glob_pattern = 'data/*/outs/'
loupe_files = sorted(glob(glob_pattern))
location_files = sorted(glob(f'{glob_pattern}/spatial/tissue_positions.csv'))

raw_filename = 'results/raw.concat.h5ad'
if not os.path.exists(raw_filename):
    adatas = list()
    spatial_uns_dict = dict()
    # spatial_uns_dict['spatial'] = dict()
    for path, loc_path in tqdm(zip(loupe_files,location_files)):

        a_ = sq.read.visium(path,  counts_file = 'filtered_feature_bc_matrix.h5') # load_images = False,
        
        df = pd.read_csv(loc_path)
        df = df.set_index('barcode')
        a_.obs = df.loc[a_.obs.index]

        a_.obs['batch'] = path
        library_id = list(a_.uns['spatial'].keys())[0]
        a_.obs['library_id'] = library_id
        a_.var_names_make_unique()

        adatas.append(a_)
        spatial_uns_dict[library_id] = a_.uns['spatial'][library_id]

    adata = anndata.concat(adatas, merge = 'same', uns_merge = 'same')
    adata.obsm['spatial'] = np.array(adata.obs[['pxl_col_in_fullres', 'pxl_row_in_fullres']])
    adata.uns['spatial'] = spatial_uns_dict

    adata.raw = adata.copy()
    adata.write(raw_filename)


adata_raw = sc.read(raw_filename)
sthd.qc.calculate_qc_metrics(adata_raw)
sthd.pl.qc_overview(adata_raw, save_name="qc_pre_filter.pdf")
sthd.io.save_h5ad(adata_raw, "output/raw_data.h5ad")

# adata.uns['spatial'] = spatial_dict
adata = adata_raw
os.makedirs('figures/qc/spatial/', exist_ok = True)
for library_id in tqdm(adata.obs['library_id'].unique()):
    a_ = adata[adata.obs['library_id'] == library_id].copy()
    sc.pl.spatial(a_, color = 'log1p_total_counts', library_id = library_id, frameon = False, show = False, title = f'{library_id}\nlog total counts') #, crop_coord = (0,60000,0,60000))
    plt.savefig(f'figures/qc/spatial/log_total_counts_{library_id}.png', bbox_inches = 'tight')
    plt.close()


annotation_files = glob('raw_data/UpdatesCSVfileSingleCell05012024/*.csv')

f2lib_id = {
    'Pathologist Annotations_SIngleCell11.csv': 'S21-10641-A4',
    'Pathologist Annotations_SingleCell10.csv': 'S21-9843-A7',
    'Pathologist Annotations_SingleCell8.csv': 'S21-5251-G3',
    'Pathologist Annotations_SingleCell17.csv': 'S21-22407-D2',
    'Pathologist Annotations_SingleCell16.csv': 'S21-21781-C15',
    'Pathologist Annotations_SingleCell2.csv': 'S21-435-A6',
    'Pathologists Annotations_SingleCell18.csv': 'S21-24036-A5',
    'Pathologist Annotations_SingleCell3.csv': 'S21-1305-A1',
}

# --- Parse and aggregate annotations ---
annotations = []
for path in annotation_files:
    fname = Path(path).name
    # Extract substring inside parentheses

    if fname not in f2lib_id.keys():
        print(f"⚠️ Could not parse library ID from: {fname}")
        continue

    # Normalize library_id: "S21 435 A6" -> "S21-435-A6"
    lib_id = f2lib_id[fname]

    # Load and standardize the annotation CSV
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    df["library_id"] = lib_id
    df.rename(columns={"Barcode": "barcode", "Pathologist Annotations": "pathologist_annotation"}, inplace=True)
    annotations.append(df)

# Concatenate all annotations into one DataFrame
anno_df = pd.concat(annotations, ignore_index=True)

# --- Normalize IDs for matching ---
anno_df["library_id"] = anno_df["library_id"].astype(str)
anno_df["barcode"] = anno_df["barcode"].astype(str)
adata.obs["library_id"] = adata.obs["library_id"].astype(str)

# --- Merge annotations into adata.obs ---
obs = adata.obs.copy()
obs["barcode"] = obs.index

merged = obs.merge(
    anno_df[["barcode", "library_id", "pathologist_annotation"]],
    how="left",
    on=["barcode", "library_id"],
)

# Assign back and restore index
adata.obs = merged.set_index("barcode")
adata.obs_names_make_unique()

# --- Verify ---
print("✅ Annotation merged successfully.")
print(adata.obs["pathologist_annotation"].value_counts(dropna=False)) #.head())



# ---- canonical cleanup ----
cleanup_dict = {
    "Non-Solid": "Non-Solid",
    "Non-Solid Blood Vessel": "Non-Solid Blood Vessel",
    "Non-Solid Bronchus": "Non-Solid Bronchus",
    "Non-Solid TLS": "Non-Solid TLS",
    "Non-solid": "Non-Solid",
    "Non-Solid Blood Vessel ": "Non-Solid Blood Vessel",

    "Normal": "Normal",
    "Normal Blood Vessel": "Normal Blood Vessel",
    "Normal Bronchus": "Normal Bronchus",
    "Normal- Bronchus": "Normal Bronchus",

    "Sold Blood Vessel": "Solid Blood Vessel",  # typo
    "Solid": "Solid",
    "Solid Blood Vessel": "Solid Blood Vessel",
    "Solid Bronchus": "Solid Bronchus",
    "Solid Scar Tissue": "Solid Scar Tissue",
    "Solid-TLS": "Solid TLS",

    "TLS Non-solid": "TLS Non-Solid",
    "TLS Normal": "TLS Normal",
    "TLS Solid": "TLS Solid",
}

adata.obs["pathologist_annotation"] = (
    adata.obs["pathologist_annotation"]
    .astype(str)
    .str.strip()
    .replace("nan", np.nan)
    .replace(cleanup_dict)
)

# ---- decomposition rules ----
def decompose_annotation(label):
    if pd.isna(label):
        return np.nan, np.nan
    label = str(label)
    # Solidity
    if "Solid" in label and "Non-Solid" not in label and "Normal" not in label:
        solidity = "Solid"
    elif "Non-Solid" in label:
        solidity = "Non-Solid"
    elif "Normal" in label:
        solidity = "Normal"
    else:
        solidity = "Unknown"

    # Architecture
    if "Vessel" in label:
        arch = "Blood Vessel"
    elif "Bronchus" in label:
        arch = "Bronchus"
    elif "TLS" in label:
        arch = "TLS"
    elif "Scar" in label:
        arch = "Scar Tissue"
    else:
        arch = "Core"  # generic tissue core
    return solidity, arch

# Apply decomposition
solidity_arch = adata.obs["pathologist_annotation"].apply(lambda x: decompose_annotation(x))
adata.obs[["solidity_type", "architecture_type"]] = pd.DataFrame(solidity_arch.tolist(), index=adata.obs.index)

# Cast to category
adata.obs["pathologist_annotation"] = adata.obs["pathologist_annotation"].astype("category")
adata.obs["solidity_type"] = adata.obs["solidity_type"].astype("category")
adata.obs["architecture_type"] = adata.obs["architecture_type"].astype("category")

print("✅ Categories — solidity:", adata.obs["solidity_type"].cat.categories.tolist())
print("✅ Categories — architecture:", adata.obs["architecture_type"].cat.categories.tolist())


sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color= 'pathologist_annotation')
sc.pl.umap(adata, color= 'pathologist_annotation', show = False)
plt.savefig('figures/umap_pathologist_annotation', bbox_inches = 'tight')
plt.close()
# --- (Optional) Save updated AnnData ---
adata.write("results/adata.annotation.h5ad")

