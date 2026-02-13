#!/usr/bin/env python3
"""
Step 1 of a two-step inventory process: extract features from listing + notebooks.

Reads the remote listing (metadata/remote_csv_listing.txt), scans all ipynb
files for references to each listed file, extracts context snippets and
heuristic columns, and writes metadata/remote_file_inventory.csv.

Step 2 (separate): A human or AI annotator reviews the CSV and the notebook
context (detailed_annotation, notebook_refs) and fills the detailed_descriptor
column in natural language. See scripts/README.md "Step 2b: Annotate
detailed_descriptor".

Output columns:
  relative_path, remote_full_path, size_bytes, parent_dir, filename,
  descriptor (short heuristic label), data_type, analysis_stage,
  detailed_annotation (raw snippets from notebooks; input for Step 2),
  notebook_refs, needs_download,
  detailed_descriptor (empty by script; fill in Step 2 with natural-language summary).

Usage (from project root):
  python scripts/build_remote_inventory.py
  python scripts/build_remote_inventory.py --listing metadata/remote_csv_listing.txt
  python scripts/build_remote_inventory.py --skip-notebooks  # fast; no context
  python scripts/build_remote_inventory.py --merge-descriptors metadata/remote_file_inventory.csv  # preserve existing detailed_descriptor when re-running
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

REMOTE_ROOT = "/home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2"

# Path prefixes in notebooks that we strip to get relative path
NOTEBOOK_PATH_PREFIXES = [
    "/athena/elementolab/scratch/dym2001/notebooks/hyperion/DLBCLv2/",
    "/home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2/",
    "/home/dym2001/athena/data/hyperion/DLBCL/",
    "DLBCLv2/",
]

# data_type inference
NORMALIZED_KEYWORDS = {"normalized", "lognorm", "log_norm", "scaled", "norm"}
PROCESSED_KEYWORDS = {
    "merged", "quant", "seurat", "cluster", "clustering", "expression",
    "counts", "stroma1_quant", "stroma2_quant", "tcell1_quant", "tcell2_quant",
    "abundance", "community",
}
RAW_KEYWORDS = {"raw", "full", "mask", "tiff", "txt"}
METADATA_KEYWORDS = {"clinical", "cellid", "cell_id", "meta", "mutation", "bcca", "seurat_cellid"}

# analysis_stage from parent dir
STAGE_PROCESSING = {"stroma_1_preprocessing", "stroma_2_preprocessing", "tcell_1_preprocessing", "tcell_2_preprocessing"}
STAGE_DOWNSTREAM = {"stroma_spatial", "total_tumor", "stroma_merged", "tcell_merged", "vessel", "tcell_myeloid", "tcell_spatial"}
STAGE_CLINICAL = {"clinical_analysis", "mutations"}


def project_root() -> Path:
    p = Path(__file__).resolve().parent
    if p.name == "scripts" and (p.parent / "Mission.md").exists():
        return p.parent
    for parent in p.parents:
        if (parent / "projects" / "imc" / "lymph_dlbcl" / "Mission.md").exists():
            return parent / "projects" / "imc" / "lymph_dlbcl"
    return p.parent


def load_remote_listing(listing_path: Path) -> list[tuple[str, int | None]]:
    """Parse remote listing. Each line: relative_path [size_bytes]. Any extension."""
    rows: list[tuple[str, int | None]] = []
    with open(listing_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split(maxsplit=1)
            rel = parts[0].strip()
            if not rel:
                continue
            size: int | None = None
            if len(parts) >= 2 and parts[1].strip().isdigit():
                size = int(parts[1].strip())
            rows.append((rel, size))
    return rows


def normalize_notebook_path(path: str) -> str:
    path = path.replace("\\", "/").strip()
    for prefix in NOTEBOOK_PATH_PREFIXES:
        if prefix in path:
            path = path.split(prefix)[-1].lstrip("/")
            break
    return path


def extract_refs_from_notebooks(
    notebooks_root: Path,
    listed_paths: set[str],
) -> dict[str, list[tuple[str, str]]]:
    """
    For each listed path, find notebook cells that reference it (by path or basename).
    Returns: relative_path -> [(notebook_name, snippet), ...]
    """
    path_to_refs: dict[str, list[tuple[str, str]]] = {p: [] for p in listed_paths}
    basename_to_paths: dict[str, list[str]] = {}
    for p in listed_paths:
        base = Path(p).name
        basename_to_paths.setdefault(base, []).append(p)

    for nb_path in notebooks_root.rglob("*.ipynb"):
        try:
            nb = json.loads(nb_path.read_text(encoding="utf-8", errors="replace"))
        except Exception:
            continue
        cells = nb.get("cells", [])
        nb_name = nb_path.name
        for cell in cells:
            source = cell.get("source", [])
            if isinstance(source, list):
                cell_lines = "".join(source).split("\n")
            else:
                cell_lines = source.split("\n")
            for line in cell_lines:
                line_stripped = line.strip()
                if not line_stripped or line_stripped.startswith("#"):
                    continue
                for ext in (".csv", ".h5ad", ".rds", ".h5"):
                    if ext not in line:
                        continue
                    for match in re.finditer(r'["\']([^"\']*' + re.escape(ext) + r')["\']', line):
                        quoted = match.group(1).replace("\\", "/").strip()
                        norm = normalize_notebook_path(quoted)
                        if not norm:
                            continue
                        snippet = line_stripped[:200] + ("..." if len(line_stripped) > 200 else "")
                        if norm in listed_paths:
                            path_to_refs[norm].append((nb_name, snippet))
                        else:
                            base = Path(norm).name
                            if base in basename_to_paths:
                                for lp in basename_to_paths[base]:
                                    path_to_refs[lp].append((nb_name, snippet))
                break  # one extension match per line is enough
    return path_to_refs


def infer_data_type(relative_path: str) -> str:
    """One of: raw | normalized | processed | metadata | unknown."""
    lower = relative_path.lower()
    parts = relative_path.replace("\\", "/").split("/")
    parent = (parts[-2] if len(parts) >= 2 else "").lower()
    fname = (parts[-1] if parts else "").lower()
    combined = parent + " " + fname
    if any(k in combined for k in NORMALIZED_KEYWORDS):
        return "normalized"
    if any(k in combined for k in METADATA_KEYWORDS):
        return "metadata"
    if any(k in combined for k in PROCESSED_KEYWORDS):
        return "processed"
    if any(k in combined for k in RAW_KEYWORDS) and "merged" not in combined:
        return "raw"
    return "unknown"


def infer_analysis_stage(relative_path: str, notebook_names: list[str]) -> str:
    """One of: processing | downstream | clinical | visualization | metadata | other."""
    parts = relative_path.replace("\\", "/").split("/")
    parent = parts[-2] if len(parts) >= 2 else ""
    parent_lower = parent.lower()
    fname_lower = (parts[-1] if parts else "").lower()
    nb_names = " ".join(notebook_names).lower()
    if parent_lower in STAGE_PROCESSING or "preprocessing" in parent_lower:
        return "processing"
    if any(d in parent_lower for d in ("stroma_spatial", "total_tumor", "stroma_merged", "tcell_merged", "vessel", "tcell_myeloid", "tcell_spatial", "1.29.23", "spatial")):
        return "downstream"
    if parent_lower in STAGE_CLINICAL or "clinical" in parent_lower or "mutation" in parent_lower:
        return "clinical"
    if "clinical" in nb_names or "survival" in nb_names or "mutation" in nb_names:
        return "clinical"
    if "viz" in nb_names or "plot" in nb_names or "figure" in nb_names or "image" in nb_names:
        return "visualization"
    if "path_sheet" in fname_lower or "metadata" in fname_lower or "cellid" in fname_lower or "clinical" in fname_lower:
        return "metadata"
    if "preprocessing" in parent_lower:
        return "processing"
    return "other"


def infer_needs_download(
    relative_path: str,
    data_type: str,
    analysis_stage: str,
    detailed_annotation: str,
) -> str:
    """One of: yes | no | review. Prioritize 'yes' for inputs we need for AnnData."""
    lower = relative_path.lower()
    fname = Path(relative_path).stem.lower()
    # First: files we need for building AnnData -> yes
    if data_type == "normalized":
        return "yes"
    if data_type == "processed" and any(k in lower for k in ("quant", "merged", "expression", "counts", "cluster", "abundance")):
        return "yes"
    if data_type == "raw" and any(k in lower for k in ("quant", "expression", "counts", "full")):
        return "yes"
    if data_type == "metadata" and any(k in lower for k in ("cellid", "cell_id", "seurat_cellid", "sample", "clinical")):
        return "yes"
    # Clear outputs -> no
    if analysis_stage == "visualization" and any(ext in lower for ext in (".png", ".jpg", ".pdf", ".svg")):
        return "no"
    if re.match(r"^\d{1,2}\.\d{1,2}\.\d{2}\.", fname) and ("write" in detailed_annotation.lower() or "to_csv" in detailed_annotation.lower()):
        return "no"
    if "write.csv(" in detailed_annotation or ".to_csv(" in detailed_annotation:
        if "read_csv" not in detailed_annotation and "read.csv" not in detailed_annotation and "fread" not in detailed_annotation:
            return "no"
    return "review"


def infer_descriptor(relative_path: str) -> str:
    parts = relative_path.replace("\\", "/").split("/")
    parent = parts[-2] if len(parts) >= 2 else ""
    fname = parts[-1] if parts else ""
    if "stroma_1" in parent or "stroma1" in relative_path.lower():
        return "stromal_panel_" + (parent or fname)
    if "stroma_2" in parent or "stroma2" in relative_path.lower():
        return "stromal_panel_" + (parent or fname)
    if "tcell_1" in parent or "tcell1" in relative_path.lower():
        return "immune_panel_" + (parent or fname)
    if "tcell_2" in parent or "tcell2" in relative_path.lower():
        return "immune_panel_" + (parent or fname)
    return parent or fname


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    root = project_root()
    default_listing = root / "metadata" / "remote_csv_listing.txt"
    parser.add_argument("--listing", type=Path, default=default_listing, help="Path to remote listing file")
    parser.add_argument("--project-root", type=Path, default=root, help="Project root")
    parser.add_argument("--remote-root", type=str, default=REMOTE_ROOT, help="Remote path prefix")
    parser.add_argument("--skip-notebooks", action="store_true", help="Skip notebook scan; annotations will be minimal")
    parser.add_argument("-o", "--output", type=Path, default=None, help="Output CSV path (default: metadata/remote_file_inventory.csv)")
    parser.add_argument(
        "--merge-descriptors",
        type=Path,
        default=None,
        help="Path to a previous inventory CSV; copy detailed_descriptor for matching relative_path so re-runs preserve annotations",
    )
    args = parser.parse_args()
    root = args.project_root
    listing_path = args.listing if args.listing.is_absolute() else root / args.listing
    metadata_dir = root / "metadata"
    notebooks_root = root / "notebooks" / "dlbcl_notebooks" / "DLBCLv2"
    out_path = args.output or metadata_dir / "remote_file_inventory.csv"

    if not listing_path.exists():
        print(f"Listing file not found: {listing_path}", file=sys.stderr)
        return 1

    rows = load_remote_listing(listing_path)
    print(f"Loaded {len(rows)} paths from {listing_path}", file=sys.stderr)

    listed_paths = set(r[0] for r in rows)

    refs = {}
    if not args.skip_notebooks and notebooks_root.exists():
        refs = extract_refs_from_notebooks(notebooks_root, listed_paths)
        total_refs = sum(len(v) for v in refs.values())
        print(f"Found {total_refs} notebook references for {len([p for p in refs if refs[p]])} files", file=sys.stderr)
    else:
        refs = {p: [] for p in listed_paths}
        if args.skip_notebooks:
            print("Skipping notebook scan (--skip-notebooks)", file=sys.stderr)

    merged_descriptors: dict[str, str] = {}
    if args.merge_descriptors:
        merge_path = args.merge_descriptors if args.merge_descriptors.is_absolute() else root / args.merge_descriptors
        if merge_path.exists():
            with open(merge_path, encoding="utf-8", newline="") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if "relative_path" in row and "detailed_descriptor" in row and (row.get("detailed_descriptor") or "").strip():
                        merged_descriptors[row["relative_path"]] = row["detailed_descriptor"].strip()
            print(f"Merged {len(merged_descriptors)} existing detailed_descriptor values from {merge_path}", file=sys.stderr)
        else:
            print(f"Warning: --merge-descriptors file not found: {merge_path}", file=sys.stderr)

    inventory = []
    seen = set()
    for rel_path, size_bytes in rows:
        if rel_path in seen:
            continue
        seen.add(rel_path)
        parts = rel_path.replace("\\", "/").split("/")
        parent_dir = "/".join(parts[:-1]) if len(parts) > 1 else ""
        filename = parts[-1] if parts else rel_path
        full_path = f"{args.remote_root.rstrip('/')}/{rel_path.lstrip('/')}"

        descriptor = infer_descriptor(rel_path)
        data_type = infer_data_type(rel_path)
        notebook_refs_list = refs.get(rel_path, [])
        notebook_names = list(dict.fromkeys([n for n, _ in notebook_refs_list]))
        analysis_stage = infer_analysis_stage(rel_path, notebook_names)

        # Detailed annotation: concatenate unique snippets (notebook: snippet)
        snippets = []
        seen_snip = set()
        for nb_name, snip in notebook_refs_list:
            key = (nb_name, snip[:150])
            if key in seen_snip:
                continue
            seen_snip.add(key)
            snippets.append(f"{nb_name}: {snip}")
        detailed_annotation = " | ".join(snippets[:5]) if snippets else ""

        needs_download = infer_needs_download(rel_path, data_type, analysis_stage, detailed_annotation)
        notebook_refs = "; ".join(notebook_names[:10])
        detailed_descriptor = merged_descriptors.get(rel_path, "")  # filled in Step 2; empty or from --merge-descriptors

        inventory.append({
            "relative_path": rel_path,
            "remote_full_path": full_path,
            "size_bytes": size_bytes if size_bytes is not None else "",
            "parent_dir": parent_dir,
            "filename": filename,
            "descriptor": descriptor,
            "data_type": data_type,
            "analysis_stage": analysis_stage,
            "detailed_annotation": detailed_annotation,
            "notebook_refs": notebook_refs,
            "needs_download": needs_download,
            "detailed_descriptor": detailed_descriptor,
        })

    inventory.sort(key=lambda r: (r["parent_dir"], r["filename"]))
    metadata_dir.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "relative_path", "remote_full_path", "size_bytes", "parent_dir", "filename",
        "descriptor", "data_type", "analysis_stage", "detailed_annotation", "notebook_refs", "needs_download",
        "detailed_descriptor",
    ]
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(inventory)
    print(f"Wrote {len(inventory)} rows to {out_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
