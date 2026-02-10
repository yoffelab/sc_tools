# -*- coding: utf-8 -*-
"""
Single JSON mapping for VisiumHD:
- Load metadata/cluster_map_with_color.json which maps leiden -> {cluster_key, cluster_color}
- Map obs['leiden'] -> obs['cluster'] and set uns['cluster_colors'] in category order
- Plot UMAP and per-library spatial for both 'leiden' and 'cluster'
- Run rank_genes_groups and matrixplots for both keys using RdBu_r
- Export DE CSVs for both keys
"""

from pathlib import Path
import json
import re
from typing import Dict, Sequence
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# Configuration
# -----------------------------
ADATA_PATH = Path("results/scvi.leiden.h5ad")
MAP_JSON = Path("metadata/cluster_map_with_color.json")
ADATA_OUT_PATH = Path("results/scvi.leiden.phenotyped.h5ad")

FIG_DIR = Path("figures")
SPATIAL_DIR = FIG_DIR / "spatial"
META_DIR = Path("metadata")
FIG_DIR.mkdir(parents=True, exist_ok=True)
SPATIAL_DIR.mkdir(parents=True, exist_ok=True)
META_DIR.mkdir(parents=True, exist_ok=True)

LEIDEN_KEY = "leiden"
CLUSTER_KEY = "cluster"
PROCESS_KEYS: Sequence[str] = (LEIDEN_KEY, CLUSTER_KEY)

# Matrixplot display parameters
CMAP = "RdBu_r"
N_TOP_GENES_FOR_HEATMAP = 5
LFC_VMIN, LFC_VMAX = -4, 4
MIN_LFC_FOR_DISPLAY = 4

DE_METHOD = "wilcoxon"  # or "wilcoxon"
USE_RAW = False

# -----------------------------
# JSON IO
# -----------------------------
def read_mapping_json(path: Path) -> Dict:
    with path.open("r") as f:
        return json.load(f)

def _is_hex_color(s: str) -> bool:
    return bool(re.fullmatch(r"#([0-9a-fA-F]{6})", s or ""))

# -----------------------------
# Mapping and color logic
# -----------------------------
def apply_leiden_mapping(
    adata,
    mapping: Dict,
    leiden_key: str = LEIDEN_KEY,
    cluster_key: str = CLUSTER_KEY,
) -> None:
    assert leiden_key in adata.obs, f"obs lacks {leiden_key}"
    if not pd.api.types.is_categorical_dtype(adata.obs[leiden_key]):
        adata.obs[leiden_key] = adata.obs[leiden_key].astype("category")

    defaults = mapping.get("_defaults", {})
    unknown_label = defaults.get("unknown_label", "Unknown")
    unknown_color = defaults.get("unknown_color", "#cccccc")

    # Validate entries and build label->color map
    label_to_color: Dict[str, str] = {}
    missing_entries = []
    color_conflicts = []

    leiden_levels = list(adata.obs[leiden_key].cat.categories)
    for lvl in leiden_levels:
        entry = mapping.get(str(lvl))
        if entry is None:
            missing_entries.append(str(lvl))
            continue
        label = entry.get("cluster_key")
        color = entry.get("cluster_color")
        if not isinstance(label, str) or not label:
            missing_entries.append(str(lvl))
            continue
        if not _is_hex_color(color):
            # fallback to unknown color if malformed
            color = unknown_color

        # Track color consistency per label
        if label in label_to_color and label_to_color[label] != color:
            color_conflicts.append((label, label_to_color[label], color))
        else:
            label_to_color[label] = color

    if missing_entries:
        print(f"Warning: missing mapping for leiden levels {missing_entries}. They will be set to '{unknown_label}'.")

    if color_conflicts:
        print("Warning: color conflicts detected for the following labels. First color kept:")
        for lab, c_old, c_new in color_conflicts:
            print(f" - {lab}: {c_old} vs {c_new}")

    # Map obs column
    def _label_for_level(x: str) -> str:
        entry = mapping.get(str(x))
        if entry is None:
            return unknown_label
        label = entry.get("cluster_key")
        return label if isinstance(label, str) and label else unknown_label

    adata.obs[cluster_key] = adata.obs[leiden_key].astype(str).map(_label_for_level).astype("category")

    # Order categories by frequency
    ordered = (
        adata.obs[cluster_key]
        .value_counts()
        .sort_values(ascending=False)
        .index
        .tolist()
    )
    adata.obs[cluster_key] = adata.obs[cluster_key].cat.reorder_categories(ordered)

    # Build palette aligned to ordered categories
    palette = []
    for label in ordered:
        color = label_to_color.get(label, unknown_color)
        palette.append(color)
    adata.uns[f"{cluster_key}_colors"] = palette

# -----------------------------
# Plotting
# -----------------------------
def plot_umap(adata, keys: Sequence[str], out_png: Path) -> None:
    sc.pl.umap(
        adata,
        color=list(keys) + ["library_id"],
        frameon=False,
        wspace=0.25,
        show=False,
    )
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

def plot_spatial_by_library(adata, color_key: str, out_dir: Path) -> None:
    assert "library_id" in adata.obs, "obs lacks library_id required for VisiumHD spatial plotting"
    for lib in adata.obs["library_id"].astype(str).unique():
        adata_sub = adata[adata.obs["library_id"].astype(str) == lib].copy()
        try:
            sc.pl.spatial(
                adata_sub,
                color=color_key,
                library_id=lib,
                frameon=False,
                show=False,
            )
            plt.savefig(out_dir / f"{color_key}_{lib}.png", dpi=300, bbox_inches="tight")
            plt.close()
        except Exception as e:
            print(f"Skip spatial plot for library {lib} and key {color_key} due to: {e}")

# -----------------------------
# Differential expression and matrixplot
# -----------------------------
def run_de(adata, groupby: str) -> None:
    assert groupby in adata.obs, f"obs lacks {groupby}"
    assert pd.api.types.is_categorical_dtype(adata.obs[groupby]), f"obs['{groupby}'] must be categorical"
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        n_genes=adata.shape[1],
        use_raw=USE_RAW,
        method=DE_METHOD,
    )

def export_de_csv(adata, groupby: str, out_csv: Path) -> None:
    cats = list(adata.obs[groupby].cat.categories)
    df = sc.get.rank_genes_groups_df(adata, group=cats)
    df.to_csv(out_csv, index=False)

def save_matrixplot(adata, groupby: str, out_png: Path) -> None:
    sc.pl.rank_genes_groups_matrixplot(
        adata,
        groupby=groupby,
        n_genes=N_TOP_GENES_FOR_HEATMAP,
        values_to_plot="logfoldchanges",
        min_logfoldchange=MIN_LFC_FOR_DISPLAY,
        vmin=LFC_VMIN,
        vmax=LFC_VMAX,
        cmap=CMAP,
        colorbar_title="log fold change",
        show=False,
    )
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

# -----------------------------
# Main
# -----------------------------
def main():
    adata = sc.read(ADATA_PATH)
    assert "library_id" in adata.obs, "obs lacks library_id"

    mapping = read_mapping_json(MAP_JSON)
    apply_leiden_mapping(adata, mapping, leiden_key=LEIDEN_KEY, cluster_key=CLUSTER_KEY)

    # UMAP and spatial for both keys
    plot_umap(adata, keys=[LEIDEN_KEY, CLUSTER_KEY], out_png=FIG_DIR / f"umap_{LEIDEN_KEY}_and_{CLUSTER_KEY}.png")
    plot_umap(adata, keys=[LEIDEN_KEY], out_png=FIG_DIR / f"umap_{LEIDEN_KEY}.png")
    plot_umap(adata, keys=[CLUSTER_KEY], out_png=FIG_DIR / f"umap_{CLUSTER_KEY}.png")
    plot_spatial_by_library(adata, color_key=LEIDEN_KEY, out_dir=SPATIAL_DIR)
    plot_spatial_by_library(adata, color_key=CLUSTER_KEY, out_dir=SPATIAL_DIR)

    # DE and matrixplots for both keys
    for key in PROCESS_KEYS:
        run_de(adata, groupby=key)
        export_de_csv(adata, groupby=key, out_csv=META_DIR / f"{key}_rank_genes.csv")
        save_matrixplot(adata, groupby=key, out_png=FIG_DIR / f"rank_gene_group_matrixplot_{key}.png")

    print("Done.")
    print(f"- UMAP: {FIG_DIR / f'umap_{LEIDEN_KEY}_and_{CLUSTER_KEY}.png'}")
    print(f"- Spatial per library in: {SPATIAL_DIR}")
    for key in PROCESS_KEYS:
        print(f"- DE CSV: {META_DIR / f'{key}_rank_genes.csv'}")
        print(f"- Matrixplot: {FIG_DIR / f'rank_gene_group_matrixplot_{key}.png'}")
    adata.write(ADATA_OUT_PATH)
if __name__ == "__main__":
    main()
