# in conda env with tia
import re
from glob import glob
from pathlib import Path

import cv2
import numpy as np
import openslide
import pandas as pd
import scanpy as sc


def parse_svs_library_id(svs_path):
    """
    Extract normalized library_id from full SVS path or filename.

    Handles examples such as:
      data/hne_images/S21-435_A6_2_121422.svs
      data/hne_images/S21_435_A6/S21_435_A6_2_121422.svs
      data/hne_images/slide_S21 435 A6/file.svs
    """
    path_str = str(svs_path)

    # Regex: capture S + digits + separator + digits + letter + digits
    match = re.search(r"(S\d{2,5})[-_\s]?(\d+)[-_]?[A-Z]\d+", path_str, re.IGNORECASE)
    if not match:
        print(f"⚠️ Could not extract library_id from: {path_str}")
        return None

    # More flexible: extract all relevant pieces
    # e.g. S21 435 A6 → S21-435-A6
    lib_match = re.search(r"S\d{2,5}[-_\s]?\d+[-_\s]?[A-Z]\d+", path_str, re.IGNORECASE)
    if not lib_match:
        print(f"⚠️ Could not fully parse library_id from: {path_str}")
        return None

    lib_id = lib_match.group(0)
    lib_id = lib_id.replace("_", "-").replace(" ", "-").upper()
    return lib_id


def build_svs_shape_dict(svs_files):
    """
    Returns {library_id: (height, width)} by reading each SVS file.
    """
    shape_dict = {}
    for svs in svs_files:
        lib_id = parse_svs_library_id(svs)
        if not lib_id:
            continue
        try:
            slide = openslide.OpenSlide(svs)
            w, h = slide.dimensions  # (width, height)
            shape_dict[lib_id] = (h, w)
            slide.close()
            print(f"📏 {lib_id}: {w}x{h}")
        except Exception as e:
            print(f"❌ Failed to open {svs}: {e}")
    return shape_dict


from math import cos, pi, sin


def hexagon_vertices(cx, cy, radius):
    """Return (x, y) vertices of a regular hexagon centered at (cx, cy)."""
    return np.array(
        [(cx + radius * cos(a), cy + radius * sin(a)) for a in np.linspace(0, 2 * pi, 7)], np.int32
    )


def draw_hex_mask(df, shape, radius, scale=1.0, label_map=None, label_col="pathologist_annotation"):
    """
    Draw a hexagonal mask with global numeric label encoding.
    label_map: {category_name: global_int_id}
    """
    mask = np.zeros(shape, dtype=np.uint16)
    r_scaled = radius * scale

    for _, row in df.iterrows():
        ann = row[label_col]
        if pd.isna(ann) or ann not in label_map:
            continue
        label_id = label_map[ann]
        cx = int(row["pxl_col_in_fullres"] * scale)
        cy = int(row["pxl_row_in_fullres"] * scale)
        hex_pts = hexagon_vertices(cx, cy, r_scaled)
        cv2.fillPoly(mask, [hex_pts], int(label_id))

    return mask


import matplotlib.cm as cm
from matplotlib.colors import ListedColormap


def save_mask_and_overlay(mask, hires_img, out_dir, lib_id, label_key, label_map):
    """
    Save mask in .npy, .png, and overlay (transparent mask over image) formats.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Define paths
    npy_path = out_dir / f"{lib_id}_{label_key}.npy"
    png_path = out_dir / f"{lib_id}_{label_key}.png"
    overlay_path = out_dir / f"{lib_id}_{label_key}_overlay.png"

    # Save as .npy
    np.save(npy_path, mask)

    # Save as .png (uint8 image)
    vmax = int(mask.max()) if mask.max() > 0 else 1
    mask_uint8 = (mask / vmax * 255).astype(np.uint8)
    cv2.imwrite(str(png_path), mask_uint8)

    # Generate overlay
    # Build color map (skip background)
    n_classes = len(label_map)
    cmap = cm.get_cmap("tab20", n_classes)
    color_mask = np.zeros_like(hires_img, dtype=np.uint8)
    for cat, idx in label_map.items():
        if idx == 0:
            continue
        color = np.array(cmap((idx - 1) / max(1, n_classes - 1))[:3]) * 255
        color_mask[mask == idx] = color.astype(np.uint8)

    # Blend with transparency
    overlay = cv2.addWeighted(hires_img, 1.0, color_mask, 0.45, 0)

    # Save overlay
    cv2.imwrite(str(overlay_path), cv2.cvtColor(overlay, cv2.COLOR_RGB2BGR))

    print(f"💾 Saved mask and overlay for {lib_id}:{label_key}")
    print(f"    ├── {npy_path.name}")
    print(f"    ├── {png_path.name}")
    print(f"    └── {overlay_path.name}")


def generate_numeric_masks(
    adata,
    svs_shape_dict,
    label_key="pathologist_annotation",
    expand_factor=1.15,
    out_dir="results/masks",
):
    """
    Generate globally consistent numeric masks for each library_id and store
    in adata.uns["masks"], and save to disk as .npy, .png, and overlay images.
    """
    adata.uns.setdefault("masks", {})

    # Global label dictionary (shared across slides)
    global_cats = list(adata.obs[label_key].cat.categories)
    label_map = {cat: i + 1 for i, cat in enumerate(global_cats)}
    print(f"✅ {label_key}: using {len(global_cats)} global categories")

    for lib in adata.obs["library_id"].unique():
        if lib not in adata.uns["spatial"]:
            print(f"⚠️ {lib} missing in adata.uns['spatial'], skipping.")
            continue

        df = adata.obs[adata.obs["library_id"] == lib].copy()
        if label_key not in df.columns:
            continue

        scales = adata.uns["spatial"][lib]["scalefactors"]
        s_hires = float(scales["tissue_hires_scalef"])
        spot_diam = float(scales["spot_diameter_fullres"])
        spot_r = spot_diam / 2.0 * expand_factor

        hires_img = adata.uns["spatial"][lib]["images"]["hires"]
        h_hires, w_hires = hires_img.shape[:2]

        mask = draw_hex_mask(
            df, (h_hires, w_hires), spot_r, s_hires, label_map, label_col=label_key
        )

        # Store in adata
        adata.uns.setdefault("masks", {}).setdefault(lib, {})[label_key] = mask

        # Save mask and overlay
        save_mask_and_overlay(mask, hires_img, out_dir, lib, label_key, label_map)

        # Report
        unique_ids = np.unique(mask)
        cats_present = [cat for cat, idx in label_map.items() if idx in unique_ids]
        print(f"📊 {lib}: mask stored ({mask.shape}), labels={cats_present}")


import matplotlib.pyplot as plt


def plot_annotation_vs_mask(
    adata,
    color_key="pathologist_annotation",
    background_color="black",
):
    """
    Plot side-by-side numeric masks vs true spatial annotations.
    Uses global color mapping from adata.uns[color_key + '_colors'].
    """
    global_cats = adata.obs[color_key].cat.categories
    global_colors = adata.uns.get(f"{color_key}_colors", None)

    if global_colors is None:
        import matplotlib.cm as cm

        global_colors = [cm.tab20(i / len(global_cats)) for i in range(len(global_cats))]
        print("⚠️ No color map found, using tab20.")
    cmap = ListedColormap([background_color] + list(global_colors))

    for lib in adata.obs["library_id"].unique():
        mask = adata.uns.get("masks", {}).get(lib, {}).get(color_key, None)
        if mask is None:
            print(f"⚠️ No mask found for {lib}.")
            continue

        adata_sub = adata[adata.obs["library_id"] == lib].copy()

        fig, axes = plt.subplots(1, 2, figsize=(12, 6))
        sc.pl.spatial(
            adata_sub,
            color=color_key,
            img_key="hires",
            library_id=lib,
            title=f"{lib} — Spatial Annotation",
            ax=axes[0],
            show=False,
            alpha=1.0,
            size=1.3,
        )
        axes[1].imshow(mask, cmap=cmap, vmin=0, vmax=len(global_cats))
        axes[1].set_title(f"{lib} — Numeric Mask")
        axes[1].axis("off")
        plt.tight_layout()
        plt.show()


def plot_all_masks(
    adata,
    mask_keys=("pathologist_annotation", "solidity_type", "architecture_type"),
    background_color="black",
):
    """
    Plot all available numeric masks (annotation, solidity, architecture) per library_id.
    - Uses consistent global colors for 'pathologist_annotation'.
    - Generates random/temporary colormaps for others.
    - Shows side-by-side comparison: spatial annotation (H&E) vs mask.
    """

    print(f"🖼️ Plotting masks for {len(mask_keys)} types: {', '.join(mask_keys)}")

    for lib in adata.obs["library_id"].unique():
        print(f"\n📊 {lib}:")

        # Collect mask keys that actually exist for this lib
        available_keys = [
            key
            for key in mask_keys
            if lib in adata.uns.get("masks", {}) and key in adata.uns["masks"][lib]
        ]
        if not available_keys:
            print(f"⚠️ No masks found for {lib}, skipping.")
            continue

        adata_sub = adata[adata.obs["library_id"] == lib].copy()
        n_keys = len(available_keys)
        fig, axes = plt.subplots(n_keys, 2, figsize=(12, 5 * n_keys))
        if n_keys == 1:
            axes = np.expand_dims(axes, 0)

        for i, key in enumerate(available_keys):
            mask = adata.uns["masks"][lib][key]
            if mask is None:
                print(f"⚠️ No mask for {lib}:{key}")
                continue

            # ---- choose color map ----
            if f"{key}_colors" in adata.uns:
                cats = list(adata.obs[key].cat.categories)
                colors = adata.uns.get(f"{key}_colors", None)
                if colors is None:
                    colors = [cm.tab20(i / len(cats)) for i in range(len(cats))]
                cmap = ListedColormap([background_color] + list(colors))
                vmax = len(cats)
            else:
                vmax = int(mask.max()) if mask.max() > 0 else 1
                cmap = cm.get_cmap("viridis", vmax + 1)

            # ---- plot scatter ----
            sc.pl.spatial(
                adata_sub,
                color=key if key in adata_sub.obs else "pathologist_annotation",
                img_key="hires",
                library_id=lib,
                title=f"{lib} — Spatial ({key})",
                ax=axes[i, 0],
                show=False,
                alpha=1.0,
                size=1.3,
            )

            # ---- plot mask ----
            axes[i, 1].imshow(mask, cmap=cmap, vmin=0, vmax=vmax)
            axes[i, 1].set_title(f"{lib} — Mask ({key})")
            axes[i, 1].axis("off")

        plt.tight_layout()
        plt.show()


# Build shape dict from real SVS files
svs_files = glob("data/hne_images/*.svs")
svs_shape_dict = build_svs_shape_dict(svs_files)

# Generate masks
adata = sc.read("results/adata.annotation.h5ad")

# Regenerate consistent numeric masks
for key in ["pathologist_annotation", "solidity_type", "architecture_type"]:
    generate_numeric_masks(adata, svs_shape_dict, label_key=key, expand_factor=1.65)

sc.pl.umap(adata, color=["pathologist_annotation", "solidity_type", "architecture_type"])

plot_all_masks(
    adata,
    mask_keys=["pathologist_annotation", "solidity_type", "architecture_type"],
    background_color="black",
)

# Regenerate consistent numeric masks and save all outputs
for key in ["pathologist_annotation", "solidity_type", "architecture_type"]:
    generate_numeric_masks(
        adata, svs_shape_dict, label_key=key, expand_factor=1.65, out_dir="results/masks"
    )

adata.write("results/adata.annotated.p2.h5ad")
