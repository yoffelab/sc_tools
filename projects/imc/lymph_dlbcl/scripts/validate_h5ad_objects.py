#!/usr/bin/env python3
"""
Phase 0.2: Validate key h5ad objects for manuscript reproduction.

Load each key Seurat-converted h5ad object, report shape, obs columns,
var names, presence of DLC_code/labels/spatial coords, and data quality.

Usage:
    python scripts/validate_h5ad_objects.py

Output:
    outputs/h5ad_validation_report.md
"""

import sys
from pathlib import Path
from datetime import datetime

import anndata as ad
import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_DIR / "results"
SEURAT_DIR = RESULTS_DIR / "seurat_converted"
OUTPUT_DIR = PROJECT_DIR / "outputs"

# Key objects to validate (relative to seurat_converted/)
KEY_OBJECTS = {
    "Immune full (T2)": "tcell_2_preprocessing/t2_SO_seurat.h5ad",
    "Stromal full (S2)": "stroma_2_preprocessing/S1_seurat_SO.h5ad",
    "Immune merged": "tcell_merged/SO2_seurat_bcell.h5ad",
    "Stromal merged": "stroma_merged/SO2_seurat_bcell.h5ad",
    "Spatial communities": "stroma_spatial/SO_k30_community_cluster.h5ad",
    "TME object": "s_tme_seurat.h5ad",
    # T2 subsets with labels
    "Immune T2 T-cell": "tcell_2_preprocessing/1.29.23.tcell_seurat.h5ad",
    "Immune T2 B-cell": "tcell_2_preprocessing/t2_seurat_bcell.h5ad",
    "Immune T2 myeloid": "tcell_2_preprocessing/t2_myeloid2_seurat.h5ad",
    "Immune T2 other": "tcell_2_preprocessing/t2_other2_seurat.h5ad",
    # S2 subsets
    "Stromal S2 B-cell": "stroma_2_preprocessing/S1_seurat_bcell.h5ad",
    "Stromal S2 other": "stroma_2_preprocessing/S1_seurat_other.h5ad",
    "Stromal S2 stroma": "stroma_2_preprocessing/S1_seurat_stroma.h5ad",
    "Stromal S2 T-cell": "stroma_2_preprocessing/S1_seurat_tcell.h5ad",
}


def validate_object(name: str, path: Path) -> dict:
    """Validate a single h5ad object and return a summary dict."""
    result = {
        "name": name,
        "path": str(path.relative_to(PROJECT_DIR)),
        "exists": path.exists(),
    }

    if not path.exists():
        result["error"] = "File not found"
        return result

    try:
        adata = ad.read_h5ad(path, backed="r")
    except Exception as e:
        result["error"] = f"Failed to load: {e}"
        return result

    result["n_obs"] = adata.n_obs
    result["n_vars"] = adata.n_vars
    result["obs_columns"] = sorted(adata.obs.columns.tolist())
    result["var_names_sample"] = adata.var_names[:10].tolist()
    result["obsm_keys"] = list(adata.obsm.keys()) if adata.obsm else []
    result["layers"] = list(adata.layers.keys()) if adata.layers else []
    result["uns_keys"] = list(adata.uns.keys()) if adata.uns else []

    # Key metadata checks
    result["has_DLC_code"] = "DLC_code" in adata.obs.columns
    result["has_labels"] = "labels" in adata.obs.columns
    result["has_meta"] = "meta" in adata.obs.columns
    result["has_orig_labels"] = "orig.labels" in adata.obs.columns
    result["has_recluster"] = "recluster" in adata.obs.columns
    result["has_seurat_clusters"] = "seurat_clusters" in adata.obs.columns
    result["has_orig_ident"] = "orig.ident" in adata.obs.columns
    result["has_community"] = "community" in adata.obs.columns
    result["has_TME"] = "TME" in adata.obs.columns
    result["has_major_group"] = "major_group" in adata.obs.columns
    result["has_spatial"] = "spatial" in adata.obsm if adata.obsm else False
    result["has_X_pca"] = "X_pca" in adata.obsm if adata.obsm else False

    # Check unique values for key columns
    if result["has_DLC_code"]:
        try:
            dlc = adata.obs["DLC_code"]
            result["n_unique_DLC"] = dlc.nunique()
            result["DLC_sample"] = sorted(dlc.unique()[:5].tolist())
        except Exception:
            result["n_unique_DLC"] = "error"

    if result["has_labels"]:
        try:
            labels = adata.obs["labels"]
            result["n_unique_labels"] = labels.nunique()
            result["label_values"] = sorted(labels.unique().tolist())
        except Exception:
            result["n_unique_labels"] = "error"

    if result["has_meta"]:
        try:
            meta = adata.obs["meta"]
            result["n_unique_meta"] = meta.nunique()
            result["meta_values"] = sorted(meta.unique().tolist())
        except Exception:
            result["n_unique_meta"] = "error"

    if result["has_community"]:
        try:
            comm = adata.obs["community"]
            result["n_unique_community"] = comm.nunique()
        except Exception:
            pass

    if result["has_orig_ident"]:
        try:
            oi = adata.obs["orig.ident"]
            result["n_unique_orig_ident"] = oi.nunique()
            result["orig_ident_sample"] = sorted(oi.unique()[:5].tolist())
        except Exception:
            pass

    # Data quality
    try:
        if hasattr(adata.X, "toarray"):
            x_sample = adata.X[:100].toarray()
        else:
            x_sample = np.array(adata.X[:100])
        result["x_min"] = float(np.nanmin(x_sample))
        result["x_max"] = float(np.nanmax(x_sample))
        result["x_mean"] = float(np.nanmean(x_sample))
        result["has_negative"] = bool(np.any(x_sample < 0))
        result["pct_zero"] = float(np.mean(x_sample == 0) * 100)
    except Exception as e:
        result["x_error"] = str(e)

    adata.file.close()
    return result


def generate_report(results: list[dict]) -> str:
    """Generate markdown validation report."""
    lines = [
        "# H5AD Validation Report",
        f"\nGenerated: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"\nValidated {len(results)} objects from `results/seurat_converted/`\n",
        "---\n",
    ]

    # Summary table
    lines.append("## Summary\n")
    lines.append("| Object | n_obs | n_vars | DLC_code | labels | meta | spatial | community |")
    lines.append("|--------|-------|--------|----------|--------|------|---------|-----------|")

    for r in results:
        if "error" in r:
            lines.append(f"| {r['name']} | ERROR | - | - | - | - | - | - |")
            continue
        lines.append(
            f"| {r['name']} | {r.get('n_obs', '?'):,} | {r.get('n_vars', '?')} "
            f"| {'Y' if r.get('has_DLC_code') else 'N'} "
            f"| {'Y' if r.get('has_labels') else 'N'} "
            f"| {'Y' if r.get('has_meta') else 'N'} "
            f"| {'Y' if r.get('has_spatial') else 'N'} "
            f"| {'Y' if r.get('has_community') else 'N'} |"
        )

    lines.append("\n---\n")

    # Detailed per-object
    lines.append("## Detailed Validation\n")

    for r in results:
        lines.append(f"### {r['name']}\n")
        lines.append(f"**Path:** `{r['path']}`\n")

        if "error" in r:
            lines.append(f"**ERROR:** {r['error']}\n")
            continue

        lines.append(f"- **Shape:** {r['n_obs']:,} obs x {r['n_vars']} vars")
        lines.append(f"- **obs columns:** {', '.join(r['obs_columns'])}")
        lines.append(f"- **var names (first 10):** {', '.join(r.get('var_names_sample', []))}")
        lines.append(f"- **obsm keys:** {', '.join(r.get('obsm_keys', [])) or 'none'}")
        lines.append(f"- **layers:** {', '.join(r.get('layers', [])) or 'none'}")

        # Key flags
        flags = []
        for key in ["has_DLC_code", "has_labels", "has_meta", "has_orig_labels",
                     "has_recluster", "has_seurat_clusters", "has_orig_ident",
                     "has_community", "has_TME", "has_major_group", "has_spatial", "has_X_pca"]:
            if r.get(key):
                flags.append(key.replace("has_", ""))
        lines.append(f"- **Key columns present:** {', '.join(flags) or 'none'}")

        # DLC info
        if r.get("has_DLC_code"):
            lines.append(f"- **Unique DLC codes:** {r.get('n_unique_DLC', '?')}")
            if "DLC_sample" in r:
                lines.append(f"  - Sample: {r['DLC_sample']}")

        # Labels
        if r.get("has_labels"):
            lines.append(f"- **Unique labels:** {r.get('n_unique_labels', '?')}")
            if "label_values" in r:
                lines.append(f"  - Values: {r['label_values']}")

        # Meta
        if r.get("has_meta"):
            lines.append(f"- **Unique meta:** {r.get('n_unique_meta', '?')}")
            if "meta_values" in r:
                lines.append(f"  - Values: {r['meta_values']}")

        # orig.ident
        if r.get("has_orig_ident"):
            lines.append(f"- **Unique orig.ident:** {r.get('n_unique_orig_ident', '?')}")
            if "orig_ident_sample" in r:
                lines.append(f"  - Sample: {r['orig_ident_sample']}")

        # Community
        if r.get("has_community"):
            lines.append(f"- **Unique communities:** {r.get('n_unique_community', '?')}")

        # Data quality
        if "x_min" in r:
            lines.append(f"- **X range:** [{r['x_min']:.3f}, {r['x_max']:.3f}], mean={r['x_mean']:.3f}")
            lines.append(f"- **Has negative values:** {r['has_negative']}")
            lines.append(f"- **% zeros (sample):** {r['pct_zero']:.1f}%")

        lines.append("")

    # Cell type label coverage
    lines.append("---\n")
    lines.append("## Cell Type Label Coverage\n")
    lines.append("Objects with `labels` column (needed for cell type mapping):\n")

    for r in results:
        if r.get("has_labels") and "label_values" in r:
            lines.append(f"- **{r['name']}**: {r['n_unique_labels']} types: {r['label_values']}")

    lines.append("\nObjects with `meta` column:\n")
    for r in results:
        if r.get("has_meta") and "meta_values" in r:
            lines.append(f"- **{r['name']}**: {r['n_unique_meta']} types: {r['meta_values']}")

    # DLC_code coverage
    lines.append("\n---\n")
    lines.append("## Sample Code (DLC_code) Coverage\n")
    for r in results:
        if r.get("has_DLC_code"):
            lines.append(f"- **{r['name']}**: {r.get('n_unique_DLC', '?')} unique samples")

    lines.append("\nObjects WITHOUT DLC_code (need mapping via orig.ident or cell ID CSVs):\n")
    for r in results:
        if not r.get("has_DLC_code") and "error" not in r:
            lines.append(f"- **{r['name']}** (orig.ident: {'Y' if r.get('has_orig_ident') else 'N'})")

    return "\n".join(lines)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Phase 0.2: Validating key h5ad objects...")
    print(f"Source: {SEURAT_DIR}")
    print()

    results = []
    for name, rel_path in KEY_OBJECTS.items():
        path = SEURAT_DIR / rel_path
        print(f"  Validating: {name} ({rel_path})...")
        result = validate_object(name, path)
        results.append(result)
        if "error" in result:
            print(f"    ERROR: {result['error']}")
        else:
            print(f"    OK: {result['n_obs']:,} obs x {result['n_vars']} vars")

    # Generate report
    report = generate_report(results)
    report_path = OUTPUT_DIR / "h5ad_validation_report.md"
    report_path.write_text(report)
    print(f"\nReport written to: {report_path}")

    # Also save as CSV for programmatic use
    csv_path = OUTPUT_DIR / "h5ad_validation_summary.csv"
    summary_rows = []
    for r in results:
        summary_rows.append({
            "name": r["name"],
            "path": r.get("path", ""),
            "n_obs": r.get("n_obs", ""),
            "n_vars": r.get("n_vars", ""),
            "has_DLC_code": r.get("has_DLC_code", False),
            "has_labels": r.get("has_labels", False),
            "has_meta": r.get("has_meta", False),
            "has_spatial": r.get("has_spatial", False),
            "has_community": r.get("has_community", False),
            "has_orig_ident": r.get("has_orig_ident", False),
            "n_unique_DLC": r.get("n_unique_DLC", ""),
            "n_unique_labels": r.get("n_unique_labels", ""),
            "error": r.get("error", ""),
        })
    pd.DataFrame(summary_rows).to_csv(csv_path, index=False)
    print(f"Summary CSV: {csv_path}")


if __name__ == "__main__":
    main()
