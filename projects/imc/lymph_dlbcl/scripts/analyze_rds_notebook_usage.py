#!/usr/bin/env python3
"""
Analyze notebook content for RDS file usage and write metadata summaries.

Scans all ipynb under notebooks/dlbcl_notebooks/DLBCLv2 for readRDS/saveRDS
and path references to .rds files. Produces:
  - metadata/rds_notebook_usage.csv: each RDS, which notebooks read/write it, inferred content, summary
  - metadata/notebook_rds_summary.csv: each notebook, which RDS it uses, brief role

Usage (from project root):
  python scripts/analyze_rds_notebook_usage.py
"""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path
from collections import defaultdict

NOTEBOOK_PATH_PREFIXES = [
    "/athena/elementolab/scratch/dym2001/notebooks/hyperion/DLBCLv2/",
    "/home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2/",
    "DLBCLv2/",
]


def project_root() -> Path:
    p = Path(__file__).resolve().parent
    if p.name == "scripts" and (p.parent / "Mission.md").exists():
        return p.parent
    for parent in p.parents:
        if (parent / "projects" / "imc" / "lymph_dlbcl" / "Mission.md").exists():
            return parent / "projects" / "imc" / "lymph_dlbcl"
    return p.parent


def normalize_rds_path(path: str) -> str:
    path = path.replace("\\", "/").strip().strip('"\'')
    for prefix in NOTEBOOK_PATH_PREFIXES:
        if prefix in path:
            path = path.split(prefix)[-1].lstrip("/")
            break
    return path


def infer_rds_content(relative_path: str) -> str:
    """Infer short description of RDS content from path/filename."""
    lower = relative_path.lower()
    parts = relative_path.replace("\\", "/").split("/")
    parent = parts[-2] if len(parts) >= 2 else ""
    fname = (parts[-1] if parts else "").replace(".rds", "")
    if "stroma_1" in parent or (parent == "" and "s_tme" in fname):
        panel = "Stromal panel (S1)"
    elif "stroma_2" in parent:
        panel = "Stromal panel (S2)"
    elif "tcell_1" in parent or "t1_" in fname:
        panel = "Immune panel (T1)"
    elif "tcell_2" in parent or "t2_" in fname:
        panel = "Immune panel (T2)"
    elif "stroma_merged" in parent:
        panel = "Stromal merged"
    elif "stroma_spatial" in parent:
        panel = "Stromal spatial"
    elif "tcell_merged" in parent:
        panel = "Immune merged"
    else:
        panel = "Other"
    if "seurat_so" in fname or "_so.rds" in lower or "so_seurat" in fname or "t1_seurat_so" in fname or "t2_so_seurat" in fname:
        return f"{panel}: Full Seurat object"
    if "bcell" in fname or "bcell" in lower:
        return f"{panel}: Seurat object (B-cell subset)"
    if "tcell" in fname or "tcell" in lower:
        return f"{panel}: Seurat object (T-cell subset)"
    if "myeloid" in fname or "myeloid" in lower:
        return f"{panel}: Seurat object (myeloid subset)"
    if "stroma" in fname and "seurat" in fname:
        return f"{panel}: Seurat object (stroma subset)"
    if "other" in fname:
        return f"{panel}: Seurat object (other subset)"
    if "k30_community" in fname or "community_cluster" in fname:
        return f"{panel}: Seurat object (k=30 community clustering)"
    if "s_tme_seurat" in fname:
        return "TME (tumor microenvironment) merged Seurat object"
    return f"{panel}: Seurat or R object"


def scan_notebooks_for_rds(notebooks_root: Path, rds_paths: set[str]) -> tuple[dict, dict]:
    """
    Returns:
      rds_to_reads: relative_path -> [(notebook_name, snippet)]
      rds_to_writes: relative_path -> [(notebook_name, snippet)]
    """
    read_pat = re.compile(r"readRDS\s*\(\s*(?:file\s*=\s*)?['\"]([^'\"]+\.rds)['\"]", re.I)
    save_pat = re.compile(r"saveRDS\s*\([^,]+,\s*(?:file\s*=\s*)?['\"]([^'\"]+\.rds)['\"]", re.I)
    rds_to_reads: dict[str, list[tuple[str, str]]] = defaultdict(list)
    rds_to_writes: dict[str, list[tuple[str, str]]] = defaultdict(list)

    for nb_path in notebooks_root.rglob("*.ipynb"):
        try:
            nb = json.loads(nb_path.read_text(encoding="utf-8", errors="replace"))
        except Exception:
            continue
        nb_name = nb_path.name
        for cell in nb.get("cells", []):
            source = cell.get("source", [])
            text = "".join(source) if isinstance(source, list) else source
            for match in read_pat.finditer(text):
                raw = match.group(1)
                norm = normalize_rds_path(raw)
                key = norm if norm in rds_paths else next((p for p in rds_paths if Path(p).name == Path(norm).name), None)
                if key:
                    snippet = match.group(0)[:120] + ("..." if len(match.group(0)) > 120 else "")
                    rds_to_reads[key].append((nb_name, snippet))
            for match in save_pat.finditer(text):
                raw = match.group(1)
                norm = normalize_rds_path(raw)
                key = norm if norm in rds_paths else next((p for p in rds_paths if Path(p).name == Path(norm).name), None)
                if key:
                    snippet = match.group(0)[:120] + ("..." if len(match.group(0)) > 120 else "")
                    rds_to_writes[key].append((nb_name, snippet))
    return dict(rds_to_reads), dict(rds_to_writes)


def main() -> int:
    root = project_root()
    notebooks_root = root / "notebooks" / "dlbcl_notebooks" / "DLBCLv2"
    rds_list_path = root / "metadata" / "files_to_download_rds.csv"
    metadata_dir = root / "metadata"

    if not rds_list_path.exists():
        print(f"RDS list not found: {rds_list_path}", file=sys.stderr)
        return 1

    rds_paths = set()
    with open(rds_list_path, encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            p = (row.get("relative_path") or "").strip()
            if p and p.endswith(".rds"):
                rds_paths.add(p)

    print(f"Loaded {len(rds_paths)} RDS paths from {rds_list_path}", file=sys.stderr)

    if not notebooks_root.exists():
        print(f"Notebooks root not found: {notebooks_root}", file=sys.stderr)
        return 1

    rds_to_reads, rds_to_writes = scan_notebooks_for_rds(notebooks_root, rds_paths)
    print(f"Found read refs for {len(rds_to_reads)} RDS; write refs for {len(rds_to_writes)} RDS", file=sys.stderr)

    # RDS-centric: metadata/rds_notebook_usage.csv
    rows = []
    for rel in sorted(rds_paths):
        reads = rds_to_reads.get(rel, [])
        writes = rds_to_writes.get(rel, [])
        read_notebooks = "; ".join(sorted(set(n for n, _ in reads))[:15])
        write_notebooks = "; ".join(sorted(set(n for n, _ in writes))[:15])
        read_snippets = " | ".join([f"{n}: {s}" for n, s in reads[:3]])
        write_snippets = " | ".join([f"{n}: {s}" for n, s in writes[:3]])
        inferred = infer_rds_content(rel)
        summary = f"Read by: {read_notebooks}. Written by: {write_notebooks}." if (reads or writes) else "No refs in scanned notebooks."
        rows.append({
            "relative_path": rel,
            "inferred_content": inferred,
            "read_by_notebooks": read_notebooks,
            "written_by_notebooks": write_notebooks,
            "read_snippets": read_snippets[:500] if read_snippets else "",
            "write_snippets": write_snippets[:500] if write_snippets else "",
            "summary": summary,
        })

    out_rds = metadata_dir / "rds_notebook_usage.csv"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    with open(out_rds, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["relative_path", "inferred_content", "read_by_notebooks", "written_by_notebooks", "read_snippets", "write_snippets", "summary"])
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} rows to {out_rds}", file=sys.stderr)

    # Notebook-centric: metadata/notebook_rds_summary.csv
    nb_to_rds: dict[str, list[tuple[str, str]]] = defaultdict(list)  # notebook -> [(rds_path, role)]
    for rel, reads in rds_to_reads.items():
        for nb_name, _ in reads:
            nb_to_rds[nb_name].append((rel, "read"))
    for rel, writes in rds_to_writes.items():
        for nb_name, _ in writes:
            nb_to_rds[nb_name].append((rel, "write"))
    nb_rows = []
    for nb_name in sorted(nb_to_rds.keys()):
        used = nb_to_rds[nb_name]
        rds_read = "; ".join(sorted(set(r for r, role in used if role == "read")))
        rds_write = "; ".join(sorted(set(r for r, role in used if role == "write")))
        nb_rows.append({
            "notebook": nb_name,
            "rds_read": rds_read,
            "rds_written": rds_write,
            "rds_count": len(set(r for r, _ in used)),
        })

    out_nb = metadata_dir / "notebook_rds_summary.csv"
    with open(out_nb, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["notebook", "rds_read", "rds_written", "rds_count"])
        w.writeheader()
        w.writerows(nb_rows)
    print(f"Wrote {len(nb_rows)} rows to {out_nb}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
