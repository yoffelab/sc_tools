"""
Signature heatmaps with versioned output (datetime-stamped PDF/PNG).

Uses sc_tools.pl for generic heatmap/clustermap logic and sc_tools.pl.save
for versioned figure paths. Run from project dir (e.g. projects/visium/ggo_visium)
so that results/ and figures/ resolve correctly.

Output: figures/manuscript/signature_heatmaps/pdf/YYDDMM.hh.mm.<basename>.pdf
        figures/manuscript/signature_heatmaps/png/YYDDMM.hh.mm.<basename>.png
"""

from pathlib import Path
import sys
from datetime import datetime

# Ensure repo root is on path when script is run from project (e.g. ggo_visium)
_script_dir = Path(__file__).resolve().parent
_repo_root = _script_dir.parent if _script_dir.name == "scripts" else _script_dir
for _ in range(5):
    if (_repo_root / "sc_tools").is_dir():
        if str(_repo_root) not in sys.path:
            sys.path.insert(0, str(_repo_root))
        break
    _repo_root = _repo_root.parent

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

from sc_tools.pl import heatmaps, save as pl_save
from sc_tools.utils import signatures as sig_utils

# -----------------------------------------------------------------------------
# Paths (relative to cwd = PROJECT)
# -----------------------------------------------------------------------------
ADATA_PATH = Path("results/adata.normalized.scored.p35.h5ad")
OUTPUT_DIR = Path("figures/manuscript/signature_heatmaps")

# -----------------------------------------------------------------------------
# Project-specific: map pathologist_annotation -> tumor_type (ggo_visium)
# -----------------------------------------------------------------------------
PATHOLOGIST_TO_TUMOR_TYPE = {
    "Solid": "Solid Tumor",
    "Non-Solid": "Non-Solid Tumor",
    "Normal": "Normal Alveolar Cells",
    "Solid Blood Vessel": "Solid Tumor",
    "Solid Bronchus": "Solid Tumor",
    "Solid Scar Tissue": "Solid Tumor",
    "Non-Solid Blood Vessel": "Non-Solid Tumor",
    "Non-Solid Bronchus": "Non-Solid Tumor",
    "Normal Blood Vessel": "Normal Alveolar Cells",
    "Normal Bronchus": "Normal Alveolar Cells",
    "Solid TLS": "Solid Tumor",
    "TLS Solid": "Solid Tumor",
    "Non-Solid TLS": "Non-Solid Tumor",
    "TLS Non-Solid": "Non-Solid Tumor",
    "TLS Normal": "Normal Alveolar Cells",
}
VALID_TUMOR_TYPES = ["Normal Alveolar Cells", "Non-Solid Tumor", "Solid Tumor"]
SOLIDITY_ORDER = ["Normal Alveolar Cells", "Non-Solid Tumor", "Solid Tumor"]
SOLIDITY_COLORS_HEX = {
    "Normal Alveolar Cells": "#66c2a5",
    "Non-Solid Tumor": "#fc8d62",
    "Solid Tumor": "#8da0cb",
}

LIBRARY_ID_COL = "library_id"
TUMOR_TYPE_COL = "tumor_type"


def _prepare_tumor_type(adata):
    """Filter to valid tumor types; create tumor_type from pathologist_annotation if needed."""
    if TUMOR_TYPE_COL not in adata.obs.columns:
        if "pathologist_annotation" not in adata.obs.columns:
            raise KeyError(
                f"Need '{TUMOR_TYPE_COL}' or 'pathologist_annotation' in adata.obs"
            )
        adata.obs[TUMOR_TYPE_COL] = pd.Categorical(
            adata.obs["pathologist_annotation"].astype(str).replace(PATHOLOGIST_TO_TUMOR_TYPE),
            categories=SOLIDITY_ORDER,
        )
    sub = adata[adata.obs[TUMOR_TYPE_COL].isin(VALID_TUMOR_TYPES)].copy()
    if isinstance(sub.obs[TUMOR_TYPE_COL].dtype, pd.CategoricalDtype):
        sub.obs[TUMOR_TYPE_COL] = sub.obs[TUMOR_TYPE_COL].cat.remove_unused_categories()
    return sub


def main():
    dt = datetime.now()
    print("=" * 60)
    print("Signature heatmaps (sc_tools, versioned output)")
    print("=" * 60)

    if not ADATA_PATH.exists():
        raise FileNotFoundError(f"AnnData not found: {ADATA_PATH}")

    adata = sc.read_h5ad(ADATA_PATH)
    adata_sub = _prepare_tumor_type(adata)
    print(f"Spots after filtering: {adata_sub.n_obs}")

    sig_columns = sig_utils.get_signature_columns(adata_sub)
    if not sig_columns:
        raise ValueError("No signature columns (sig:..._z) found. Run score_gene_signatures first.")

    if LIBRARY_ID_COL not in adata_sub.obs.columns:
        raise KeyError(f"Missing obs column: {LIBRARY_ID_COL}")

    annotation_cols = {"Patient": LIBRARY_ID_COL, "Solidity": TUMOR_TYPE_COL}
    category_orders = {"Solidity": SOLIDITY_ORDER}

    configs = [
        ("signature_heatmap_patient_solidity", ["Patient", "Solidity"], False),
        ("signature_heatmap_solidity_patient", ["Solidity", "Patient"], False),
        ("signature_clustermap_patient_solidity", ["Patient", "Solidity"], True),
        ("signature_clustermap_solidity_patient", ["Solidity", "Patient"], True),
    ]

    for basename, sort_by, cluster in configs:
        print(f"  Building {basename}...")
        fig, _ = heatmaps.signature_score_heatmap(
            adata_sub,
            sig_columns,
            annotation_cols=annotation_cols,
            sort_by=sort_by,
            category_orders=category_orders,
            cluster=cluster,
            solidity_colors_hex=SOLIDITY_COLORS_HEX,
            legend_title="Solidity",
        )
        pdf_path, png_path = pl_save.save_figure(
            fig, basename, OUTPUT_DIR, dt=dt, create_pdf_png_folders=True
        )
        plt.close(fig)
        print(f"    -> {pdf_path.relative_to(OUTPUT_DIR)}")
        print(f"    -> {png_path.relative_to(OUTPUT_DIR)}")

    print("=" * 60)
    print(f"Done. Version prefix: {dt.strftime('%y%d%m.%H.%M')}")
    print(f"Output dir: {OUTPUT_DIR}/pdf/ and {OUTPUT_DIR}/png/")
    print("=" * 60)


if __name__ == "__main__":
    main()
