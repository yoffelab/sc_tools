#!/usr/bin/env python3
"""
Step 2b: Fill detailed_descriptor from inventory columns and write annotated CSV.

Reads metadata/remote_file_inventory.csv (script output; detailed_descriptor empty),
generates a natural-language summary for each row from descriptor, data_type,
analysis_stage, detailed_annotation, and path/filename, and writes
metadata/remote_file_inventory_annotated.csv. Keeps the original file unchanged
for traceability.

Usage (from project root):
  python scripts/annotate_inventory_descriptors.py
  python scripts/annotate_inventory_descriptors.py -i metadata/remote_file_inventory.csv -o metadata/remote_file_inventory_annotated.csv
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


def project_root() -> Path:
    p = Path(__file__).resolve().parent
    if p.name == "scripts" and (p.parent / "Mission.md").exists():
        return p.parent
    for parent in p.parents:
        if (parent / "projects" / "imc" / "lymph_dlbcl" / "Mission.md").exists():
            return parent / "projects" / "imc" / "lymph_dlbcl"
    return p.parent


def describe_one(row: dict) -> str:
    """Generate a short natural-language detailed_descriptor from inventory row."""
    rel = row.get("relative_path", "")
    parent = row.get("parent_dir", "")
    fname = row.get("filename", rel)
    descriptor = row.get("descriptor", "")
    data_type = (row.get("data_type") or "").strip().lower()
    stage = (row.get("analysis_stage") or "").strip().lower()
    ann = (row.get("detailed_annotation") or "").lower()
    needs = (row.get("needs_download") or "").strip().lower()
    lower = rel.lower()

    # Panel and content type from path
    is_stromal = "stroma" in lower and "tcell" not in lower
    is_immune = "tcell" in lower or "immune" in lower
    panel = "stromal" if is_stromal else ("immune" if is_immune else "")

    parts = []

    # Expression / count matrices (cells as observations)
    if "quant_merged" in lower or "quant_clustering" in lower:
        p = "Stromal" if is_stromal else "Immune"
        if "stroma_1" in parent or "stroma1" in lower:
            p = "Stromal panel (S1)"
        elif "stroma_2" in parent or "stroma2" in lower:
            p = "Stromal panel (S2)"
        elif "tcell_1" in parent or "tcell1" in lower:
            p = "Immune panel (T1)"
        elif "tcell_2" in parent or "tcell2" in lower:
            p = "Immune panel (T2)"
        if data_type == "normalized":
            parts.append(f"{p} normalized expression matrix; cells as rows, markers as columns.")
        else:
            parts.append(f"{p} merged expression matrix; cells as rows, markers as columns. Processed counts from prior pipeline.")
        if "read_csv" in ann or "read.csv" in ann or "fread" in ann:
            parts.append("Used as main input for clustering, spatial analysis, and downstream notebooks.")
        return " ".join(parts)

    # Cell IDs / cell metadata
    if "cellid" in lower or "cell_id" in lower or "seurat_cellid" in lower:
        p = "Stromal" if is_stromal else "Immune"
        if "S1_" in fname or "stroma_1" in parent:
            p = "Stromal panel (S1)"
        elif "S2_" in fname or "stroma_2" in parent:
            p = "Stromal panel (S2)"
        parts.append(f"{p} cell ID table mapping cells to samples or slides.")
        parts.append("Used to join cell-level data with sample metadata in downstream analysis.")
        return " ".join(parts)

    # Cluster labels / TME clusters / cluster averages
    if "cluster" in lower and ("c8" in lower or "c10" in lower or "stroma_c" in lower or "tme" in lower):
        if "avg" in lower or "average" in lower:
            parts.append("Cluster average expression or summary (e.g. C10 stromal clusters). Written by case clustering; read by spatial_plots.")
        else:
            parts.append("TME or stromal cluster assignments per cell or per sample.")
            if "viz" in ann or "visualization" in stage:
                parts.append("Read in ID and community visualization notebooks.")
        return " ".join(parts)

    # Spatial counts (large matrices)
    if "spatial_counts" in lower:
        if "norm" in lower:
            parts.append("Stromal (S1) spatial normalized count matrix; cells as rows. Used in community clustering and spatial analysis.")
        elif "merged" in lower or "bmerged" in lower:
            parts.append("Stromal (S1) spatial merged count matrix (B-cell merged). Used in spatial and community clustering.")
        else:
            parts.append("Stromal (S1) spatial count matrix; cells as rows. Large expression/count table for spatial analysis.")
        return " ".join(parts)

    # Abundance / proportions
    if "abundance" in lower or "percent" in lower and "meta" in lower:
        parts.append("Abundance or proportion table (e.g. cell type or marker per sample).")
        if "CD206" in rel or "CD163" in rel:
            parts.append("Macrophage marker abundance. Used in case clustering and RNA integration.")
        else:
            parts.append("Used in heatmaps, case clustering, or meta-analysis.")
        return " ".join(parts)

    # Path / image sheets (metadata for locating files)
    if "path_sheet" in lower or "file_path" in lower or "image_sheet" in lower:
        parts.append("Table of file paths or image paths for stroma/panels. Used to load images or link samples to file locations in visualization and renormalization.")
        return " ".join(parts)

    # Clinical / IHC outputs
    if "clinical" in stage or "IHC" in fname or "clinical" in lower:
        if "write.csv" in ann and "read" not in ann:
            parts.append("Output table from case clustering (e.g. IHC or clinical subset). Written by DLBCL_case_clustering; read by ML_framework.")
        else:
            parts.append("Clinical or IHC-derived table. Used in survival, ML framework, or case-level analysis.")
        return " ".join(parts)

    # TME table (cell-level table used across notebooks)
    if "tme_table" in lower or "tme_filt" in lower or "s_tme" in lower:
        parts.append("TME (tumor microenvironment) cell-level table: sample/cell IDs and possibly phenotype. Used in stroma_panel_tme, case clustering, and complex heatmaps.")
        return " ".join(parts)

    # B-cell / correlation tables
    if "bcell" in lower and ("table" in lower or "percent" in lower):
        parts.append("B-cell proportion or normalized table for stromal/immune panel. Written by tcell_myeloid notebooks; read by tumor complex heatmaps.")
        return " ".join(parts)

    # Community / spatial downstream
    if "community" in lower and "path" in lower:
        parts.append("Community or pairwise analysis path table. Output of pairwise viz; read by stroma_maps for visualization.")
        return " ".join(parts)

    # Antigen / tumor antigen figures
    if "antigen" in lower or "merged_antigens" in lower:
        parts.append("Merged antigen average scores or totals. Used in tumor antigen figures and case clustering.")
        return " ".join(parts)

    # Coordinate or small output (full_coo)
    if "full_coo" in lower:
        parts.append("Coordinate or small table output from tumor_antigen_figures. Likely spatial or plot-related.")
        return " ".join(parts)

    # Immune / stroma average scores
    if "immune_avg" in lower or "avg_scores" in lower:
        parts.append("Average immune or signature scores per sample. Read by complex heatmaps and tumor complex heatmaps.")
        return " ".join(parts)

    # percent_meta
    if "percent_meta" in lower:
        parts.append("Percentage or proportion metadata. Read by complex heatmaps for visualization.")
        return " ".join(parts)

    # stroma_C10 (root-level cluster file)
    if "stroma_c10" in lower or "stroma_c8" in lower:
        parts.append("Stromal TME cluster assignments (C8/C10). Used in case clustering and Hyperion RNA integration.")
        return " ".join(parts)

    # General: output of a notebook (write.csv / to_csv)
    if "write.csv" in ann or "to_csv" in ann:
        if "read" not in ann and "fread" not in ann:
            parts.append("Output file written by a notebook (e.g. filtered subset or summary table). Typically not required for building AnnData from scratch.")
            return " ".join(parts)

    # General: from descriptor and data_type
    if data_type == "normalized":
        parts.append(f"Normalized data. {descriptor or fname}. Used in downstream analysis.")
    elif data_type == "processed":
        parts.append(f"Processed output. {descriptor or fname}. Used in clustering or downstream notebooks.")
    elif data_type == "metadata":
        parts.append(f"Metadata or mapping table. {descriptor or fname}. Used to join or annotate cell/sample data.")
    elif data_type == "raw":
        parts.append(f"Raw or unprocessed data. {descriptor or fname}. May be used as input to normalization or visualization.")
    else:
        if stage == "processing":
            parts.append(f"Processing-stage file. {descriptor or fname}. Check detailed_annotation for usage.")
        elif stage == "downstream":
            parts.append(f"Downstream analysis file. {descriptor or fname}. Check detailed_annotation for usage.")
        elif stage == "visualization":
            parts.append(f"Used in visualization or mapping. {descriptor or fname}.")
        elif stage == "clinical":
            parts.append(f"Clinical or mutation-related. {descriptor or fname}.")
        else:
            parts.append(f"{descriptor or fname}. See detailed_annotation and notebook_refs for context.")
    return " ".join(parts)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    root = project_root()
    parser.add_argument("-i", "--input", type=Path, default=root / "metadata" / "remote_file_inventory.csv", help="Input inventory CSV")
    parser.add_argument("-o", "--output", type=Path, default=root / "metadata" / "remote_file_inventory_annotated.csv", help="Output annotated CSV")
    args = parser.parse_args()
    inp = args.input if args.input.is_absolute() else root / args.input
    out = args.output if args.output.is_absolute() else root / args.output

    if not inp.exists():
        print(f"Input not found: {inp}", file=sys.stderr)
        return 1

    rows = []
    with open(inp, encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        for row in reader:
            row["detailed_descriptor"] = describe_one(row)
            rows.append(row)

    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} rows to {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
