# Project Architecture: Two-Panel IMC DLBCL (lymph_dlbcl)

This document describes the **project-specific** directory layout, data sources, and data flow for the DLBCL two-panel IMC project. Repository-wide conventions are in the root `Architecture.md`.

---

## 1. Directory Layout

```text
projects/imc/lymph_dlbcl/
├── data/              # Local copy of raw or staged data (when used)
│   └── downloaded/   # Required processed/normalized/filtered files pulled from remote (Step 2)
├── figures/           # QC, manuscript, and exploratory figures
├── metadata/          # Project metadata and remote inventory
│   ├── remote_file_inventory.csv           # Script output (Step 2a); detailed_descriptor empty
│   ├── remote_file_inventory_annotated.csv # Annotated inventory with detailed_descriptor (Step 2b); use for download decisions
│   ├── remote_csv_listing.txt     # Raw find output from remote (reproducible input)
│   ├── files_to_download_rds.csv  # List of RDS paths to download (47 files)
│   ├── rds_notebook_usage.csv     # Per-RDS: which notebooks read/write it, inferred content (analyze_rds_notebook_usage.py)
│   ├── notebook_rds_summary.csv  # Per-notebook: which RDS it reads/writes (analyze_rds_notebook_usage.py)
│   ├── files_to_download.csv     # Step 2: list of remote paths to download (required for AnnData)
│   ├── download_manifest.csv     # After download: remote_path, local_path, date (traceability)
│   ├── sample_metadata.csv        # Phase 2: sample → clinical (when ready)
│   └── celltype_map.json         # Phase 4: cluster_id → celltype (when ready)
├── notebooks/         # Existing DLBCL notebooks (reference + extraction)
│   └── dlbcl_notebooks/
│       └── DLBCLv2/   # Two-panel immune and stromal preprocessing and analyses
├── results/           # AnnData checkpoints (standard names per root Architecture.md)
├── outputs/           # Intermediate outputs (logs, temp)
├── scripts/           # Project-specific ingestion and pipeline scripts
├── Mission.md
├── Architecture.md   # This file
└── Journal.md
```

---

## 2. Data Sources and Paths

### Remote (reference; no bulk download required for planning)

| Purpose | Path |
|--------|------|
| **Remote project root** | `/home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2` |
| **Access** | SSH (e.g. `ssh cayuga`) then run `find` or list commands. |

All CSV and directory-structure discovery is done against this path. Notebooks from the prior analyst may reference other paths (e.g. `/athena/elementolab/scratch/...`, `/home/dym2001/...`); the **canonical** source for this project is the backup path above.

### Local

| Purpose | Path |
|--------|------|
| **Project root** | `projects/imc/lymph_dlbcl/` |
| **Notebooks (reference)** | `projects/imc/lymph_dlbcl/notebooks/dlbcl_notebooks/DLBCLv2/` |
| **Inventory and metadata** | `projects/imc/lymph_dlbcl/metadata/` |
| **Checkpoint AnnData** | `projects/imc/lymph_dlbcl/results/` |

---

## 3. Checkpoint Nomenclature (alignment with root Architecture.md)

Once Phase 1–4 are implemented, use the standard checkpoint names:

| Phase | Standard path (under `results/`) | Description |
|-------|----------------------------------|-------------|
| 1 | `adata.raw.p1.h5ad` | Raw AnnData after ingestion and QC. |
| 2 | `adata.annotated.p2.h5ad` | With sample/clinical metadata attached. |
| 3 | `adata.normalized.p3.h5ad` | Filtered, normalized, clustered. |
| 3.5b | `adata.normalized.scored.p35.h5ad` | Normalized + signature scores (and optional automated cell typing). |
| 4 | `adata.celltyped.p4.h5ad` | Manual cell type labels applied. |

For **two-panel** structure, keep **immune** and **stromal** panels **separate**. Use panel-prefixed checkpoints (e.g. `results/adata.immune.raw.p1.h5ad`, `results/adata.stromal.raw.p1.h5ad`, `results/adata.immune.normalized.p3.h5ad`, `results/adata.stromal.normalized.p3.h5ad`). Document any additional naming in this file and in Mission.md.

---

## 4. Data Flow (Planning and Execution)

**Priority:** Reuse existing processed data from the previous analyst when the inventory shows it is suitable. Process raw only when needed.

1. **Step 1 — Inventory:** Remote `find` for `*.csv` (and optionally `*.h5ad`, `*.rds`, `*.h5`); notebook scan for path references; build `metadata/remote_csv_inventory.csv` with descriptor, inferred_role, and processed/raw tag.
2. **Step 2 (priority) — Identify, download, reuse:** From inventory, identify **required** processed/normalized/filtered files; record in `metadata/files_to_download.csv`; **download** them to `data/downloaded/`; record in `metadata/download_manifest.csv`; build AnnData from downloaded files into `results/adata.immune.*.h5ad` and `results/adata.stromal.*.h5ad`.
3. **Step 3 (fallback) — Raw:** For panel(s) without suitable processed data, ingest raw and run Phase 1; write panel-specific `results/adata.{immune|stromal}.raw.p1.h5ad` and downstream checkpoints. Keep immune and stromal pipelines separate.

---

## 5. References

- Root `Architecture.md`: checkpoint names, required metadata, phase details.
- Root `Mission.md`: toolkit phasing and entry points.
- This project `Mission.md`: roadmap and active tasks (inventory first, then reuse processed data; raw only if needed). Panels: immune and stromal, kept separate.
