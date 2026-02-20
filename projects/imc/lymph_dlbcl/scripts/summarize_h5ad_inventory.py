#!/usr/bin/env python3
"""
Summarize converted DLBCL h5ad files: obs/var counts, key columns (DLC_code, labels, meta),
obsm/layers, and inferred role (raw vs processed, panel, usefulness).

Run from project root: projects/imc/lymph_dlbcl/
  python scripts/summarize_h5ad_inventory.py [h5ad_root]

Output: metadata/h5ad_inventory_summary.csv, metadata/h5ad_inventory_summary.md
"""

from __future__ import annotations

import argparse
import csv
import sys
import warnings
from pathlib import Path

# Suppress anndata compat warnings when opening Seurat-converted h5ad
warnings.filterwarnings("ignore", message="Moving element from .uns\\['neighbors'\\]")

try:
    import anndata
except ImportError:
    print("Error: anndata is required. Install with: pip install anndata", file=sys.stderr)
    sys.exit(1)

try:
    import h5py
except ImportError:
    h5py = None  # fallback to anndata full read


def _find_h5ad_files(root: Path) -> list[Path]:
    return sorted(root.rglob("*.h5ad"))


def _relpath(path: Path, base: Path) -> str:
    try:
        return path.relative_to(base).as_posix()
    except ValueError:
        return path.as_posix()


def _infer_stage(rel: str) -> str:
    if rel.startswith("stroma_1_preprocessing") or rel.startswith("stroma_2_preprocessing"):
        return "stroma_preprocessing"
    if rel.startswith("tcell_1_preprocessing") or rel.startswith("tcell_2_preprocessing"):
        return "tcell_preprocessing"
    if rel.startswith("stroma_merged") or rel.startswith("tcell_merged"):
        return "merged"
    if "spatial" in rel or "k30" in rel.lower():
        return "spatial"
    if rel.startswith("s_tme"):
        return "TME_merged"
    return "other"


def _infer_panel(rel: str) -> str:
    if "stroma" in rel:
        return "stromal"
    if "tcell" in rel or "t1_" in rel or "t2_" in rel:
        return "immune"
    if "s_tme" in rel or "TME" in rel:
        return "TME"
    return "unknown"


def _load_rds_metadata(project_root: Path) -> dict[str, str]:
    """Load inferred_content from rds_notebook_usage by relative_path (rds)."""
    csv_path = project_root / "metadata" / "rds_notebook_usage.csv"
    out: dict[str, str] = {}
    if not csv_path.exists():
        return out
    with open(csv_path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        if "relative_path" not in (r.fieldnames or []):
            return out
        for row in r:
            rel = (row.get("relative_path") or "").strip()
            if rel.endswith(".rds"):
                # key for matching h5ad: same path with .h5ad
                key = rel[:-4] + ".h5ad"
                out[key] = (row.get("inferred_content") or "").strip()
    return out


def _fast_scan_h5py(h5ad_path: Path) -> dict | None:
    """Use h5py to read only structure (no anndata full load). Returns dict with n_obs, n_vars, obs_columns, etc. or None."""
    if not h5py:
        return None
    try:
        with h5py.File(h5ad_path, "r") as f:
            # X is (n_obs, n_vars) in AnnData
            if "X" in f:
                x = f["X"]
                if hasattr(x, "shape") and len(x.shape) >= 2:
                    n_obs, n_vars = x.shape[0], x.shape[1]
                else:
                    n_obs = f.attrs.get("n_obs", None)
                    n_vars = f.attrs.get("n_vars", None)
            else:
                n_obs = f.attrs.get("n_obs", None)
                n_vars = f.attrs.get("n_vars", None)
            obs_cols = []
            if "obs" in f:
                obs = f["obs"]
                if hasattr(obs, "keys"):
                    # obs is a group; column names are often keys (pandas 0.23+ format)
                    obs_cols = [k for k in obs.keys() if not k.startswith("_")]
                elif hasattr(obs, "dtype") and obs.dtype.names:
                    obs_cols = list(obs.dtype.names)
            obsm_keys = list(f["obsm"].keys()) if "obsm" in f and hasattr(f["obsm"], "keys") else []
            layer_keys = list(f["layers"].keys()) if "layers" in f and hasattr(f["layers"], "keys") else []
            uns_keys = list(f["uns"].keys()) if "uns" in f and hasattr(f["uns"], "keys") else []
            return {
                "n_obs": n_obs,
                "n_vars": n_vars,
                "obs_columns": obs_cols,
                "obsm_keys": obsm_keys,
                "layer_names": layer_keys,
                "uns_keys": uns_keys,
            }
    except Exception:
        return None


def summarize_one(h5ad_path: Path, project_root: Path, rds_meta: dict[str, str], h5ad_root: Path) -> dict:
    rel = _relpath(h5ad_path, h5ad_root)
    row: dict = {
        "relative_path": rel,
        "n_obs": None,
        "n_vars": None,
        "obs_columns": "",
        "has_DLC_code": False,
        "has_labels": False,
        "has_meta": False,
        "label_like_columns": "",
        "obsm_keys": "",
        "layer_names": "",
        "uns_keys": "",
        "n_unique_DLC": None,
        "inferred_content": rds_meta.get(rel, ""),
        "stage": _infer_stage(rel),
        "panel": _infer_panel(rel),
        "usefulness_notes": "",
    }
    fast = _fast_scan_h5py(h5ad_path) if h5py else None
    if fast:
        row["n_obs"] = fast.get("n_obs")
        row["n_vars"] = fast.get("n_vars")
        obs_cols = fast.get("obs_columns") or []
        row["obs_columns"] = "|".join(sorted(obs_cols))
        row["has_DLC_code"] = "DLC_code" in obs_cols
        row["has_labels"] = "labels" in obs_cols
        row["has_meta"] = "meta" in obs_cols
        label_like = [c for c in obs_cols if "label" in c.lower() or "meta" in c.lower() or "celltype" in c.lower() or "cluster" in c.lower()]
        row["label_like_columns"] = "|".join(sorted(label_like)) if label_like else ""
        row["obsm_keys"] = "|".join(sorted(fast.get("obsm_keys") or []))
        row["layer_names"] = "|".join(sorted(fast.get("layer_names") or []))
        row["uns_keys"] = "|".join(sorted(fast.get("uns_keys") or []))
    else:
        try:
            ad = anndata.read_h5ad(h5ad_path, backed="r")
            row["n_obs"] = ad.n_obs
            row["n_vars"] = ad.n_vars
            obs_cols = list(ad.obs.columns)
            row["obs_columns"] = "|".join(sorted(obs_cols))
            row["has_DLC_code"] = "DLC_code" in ad.obs.columns
            row["has_labels"] = "labels" in ad.obs.columns
            row["has_meta"] = "meta" in ad.obs.columns
            label_like = [c for c in obs_cols if "label" in c.lower() or "meta" in c.lower() or "celltype" in c.lower() or "cluster" in c.lower()]
            row["label_like_columns"] = "|".join(sorted(label_like)) if label_like else ""
            row["obsm_keys"] = "|".join(sorted(ad.obsm.keys())) if ad.obsm else ""
            row["layer_names"] = "|".join(sorted(ad.layers.keys())) if ad.layers else ""
            row["uns_keys"] = "|".join(sorted(ad.uns.keys())) if ad.uns else ""
            if row["has_DLC_code"]:
                row["n_unique_DLC"] = ad.obs["DLC_code"].nunique()
            ad.file.close()
        except Exception as e:
            row["usefulness_notes"] = f"error: {e}"
            return row

    # Heuristic: full SO = more raw; merged + cell-type = more processed; DLC_code = sample-level; labels/meta = cell-type annotated
    notes = []
    if "preprocessing" in row["stage"] and "SO" in Path(rel).stem.upper() and "bcell" not in rel and "tcell" not in rel and "stroma" not in rel and "myeloid" not in rel and "other" not in rel:
        notes.append("full_object_likely_raw")
    if row["stage"] == "merged":
        notes.append("merged_downstream")
    if row["has_DLC_code"]:
        notes.append("has_sample_code")
    if row["has_labels"] or row["has_meta"] or row["label_like_columns"]:
        notes.append("has_celltype_or_labels")
    if row["stage"] == "spatial":
        notes.append("spatial_clustering")
    row["usefulness_notes"] = "; ".join(notes) if notes else ""
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize DLBCL h5ad inventory (obs, var, DLC_code, labels, etc.)")
    parser.add_argument(
        "h5ad_root",
        nargs="?",
        default="results/seurat_converted",
        help="Root directory to search for .h5ad files (default: results/seurat_converted)",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="metadata",
        help="Directory for output CSV and MD (default: metadata)",
    )
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parent.parent
    h5ad_root = project_root / args.h5ad_root
    output_dir = project_root / args.output_dir

    if not h5ad_root.exists():
        print(f"h5ad root does not exist: {h5ad_root}; writing empty summary.", file=sys.stderr)
        output_dir.mkdir(parents=True, exist_ok=True)
        csv_path = output_dir / "h5ad_inventory_summary.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=[
                "relative_path", "n_obs", "n_vars", "n_unique_DLC", "has_DLC_code", "has_labels", "has_meta",
                "label_like_columns", "obs_columns", "obsm_keys", "layer_names", "uns_keys",
                "inferred_content", "stage", "panel", "usefulness_notes",
            ])
            w.writeheader()
        (output_dir / "h5ad_inventory_summary.md").write_text(
            "# DLBCL h5ad inventory summary\n\nNo h5ad files found. Run batch_seurat_to_h5ad.sh first, then re-run this script.\n",
            encoding="utf-8",
        )
        print(f"Wrote empty summary to {csv_path}")
        return

    paths = _find_h5ad_files(h5ad_root)
    if not paths:
        print(f"No .h5ad files found under {h5ad_root}", file=sys.stderr)
        output_dir.mkdir(parents=True, exist_ok=True)
        csv_path = output_dir / "h5ad_inventory_summary.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=[
                "relative_path", "n_obs", "n_vars", "n_unique_DLC", "has_DLC_code", "has_labels", "has_meta",
                "label_like_columns", "obs_columns", "obsm_keys", "layer_names", "uns_keys",
                "inferred_content", "stage", "panel", "usefulness_notes",
            ])
            w.writeheader()
        print(f"Wrote empty summary to {csv_path}")
        return

    rds_meta = _load_rds_metadata(project_root)
    rows = [summarize_one(p, project_root, rds_meta, h5ad_root) for p in paths]

    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "h5ad_inventory_summary.csv"
    fieldnames = [
        "relative_path", "n_obs", "n_vars", "n_unique_DLC", "has_DLC_code", "has_labels", "has_meta",
        "label_like_columns", "obs_columns", "obsm_keys", "layer_names", "uns_keys",
        "inferred_content", "stage", "panel", "usefulness_notes",
    ]
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    md_path = output_dir / "h5ad_inventory_summary.md"
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("# DLBCL h5ad inventory summary\n\n")
        f.write("Generated by `scripts/summarize_h5ad_inventory.py`. Use this to identify ")
        f.write("which files are full (raw) vs subset/processed, which have sample codes (DLC_code), ")
        f.write("and which have cell-type labels/meta.\n\n")
        f.write("## Summary table (samples × vars, key columns)\n\n")
        f.write("| relative_path | n_obs | n_vars | n_unique_DLC | DLC_code | labels | meta | label_like | stage | panel | usefulness_notes |\n")
        f.write("| --- | ---: | ---: | ---: | --- | --- | --- | --- | --- | --- | --- |\n")
        for r in rows:
            label_preview = (r["label_like_columns"] or "")[:50].replace("|", ",")
            f.write(f"| {r['relative_path']} | {r['n_obs'] or ''} | {r['n_vars'] or ''} | {r['n_unique_DLC'] or ''} | ")
            f.write(f"{'Y' if r['has_DLC_code'] else ''} | {'Y' if r['has_labels'] else ''} | {'Y' if r['has_meta'] else ''} | ")
            f.write(f"{label_preview} | {r['stage']} | {r['panel']} | {r['usefulness_notes']} |\n")
        f.write("\n## Which files are more useful / more recent?\n\n")
        f.write("- **Full (raw) objects** (one per stromal/immune panel): Look for `stage=stroma_preprocessing` or ")
        f.write("`tcell_preprocessing` with filename containing `_SO` or `seurat_SO` (e.g. S1_seurat_SO, t1_SO_seurat, t2_SO_seurat). ")
        f.write("These have the most cells before subsetting by cell type.\n")
        f.write("- **Sample code (DLC_code)**: Filter rows with `has_DLC_code=Y` to get objects that include sample/slide identifiers.\n")
        f.write("- **Cell-type annotated**: Filter by `has_labels=Y` or `has_meta=Y` or non-empty `label_like_columns` for downstream typing.\n")
        f.write("- **Merged / downstream**: `stage=merged` or `stage=spatial` are later pipeline steps (e.g. tumor scoring, spatial clustering).\n")
        f.write("\n## Full obs columns and obsm/layers (see CSV)\n\n")
        f.write("For full `obs_columns`, `obsm_keys`, `layer_names`, and `inferred_content`, open `h5ad_inventory_summary.csv`.\n")

    print(f"Summarized {len(rows)} h5ad files.")
    print(f"CSV: {csv_path}")
    print(f"MD:  {md_path}")


if __name__ == "__main__":
    main()
