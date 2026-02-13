# Mission: Two-Panel IMC DLBCL (lymph_dlbcl)

**Project:** `projects/imc/lymph_dlbcl`  
**Current Status:** Inventory complete; next: identify required files, download, build AnnData  
**Author:** Junbum Kim  
**Last Updated:** 2025-02-12

This file holds **project-specific** goals. Repository-level pipeline and phase definitions are in the root `Mission.md` and `Architecture.md`.

---

## 1. Objective

- **Scientific:** Two-panel Hyperion IMC project on DLBCL (diffuse large B-cell lymphoma). The two panels are **immune** and **stromal**; keep them **separate** (e.g. separate h5ad per panel). Produce multiple h5ad files aligned with sc_tools IMC conventions and other IMC projects in this repo.
- **Operational (priority order):** (1) Build file inventory and **identify any already-processed data** from the previous analyst; **reuse it** if suitable. (2) **Only if** suitable processed outputs are not found, process raw data and run the sc_tools pipeline (ingestion → preprocessing → checkpoint AnnData). Do not re-process from raw when existing processed data can be picked up.

---

## 2. Context and Conventions

- **Remote data (read-only reference):** `/home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2` (cayuga or equivalent).
- **Local project root:** `./projects/imc/lymph_dlbcl/`.
- **Existing assets:** Processing notebooks are in `notebooks/dlbcl_notebooks/DLBCLv2/` (two panels: **immune**, **stromal**; preprocessing per panel, merged analyses, clinical, spatial, vessel, etc.).
- **Target:** Standard checkpoint nomenclature per root `Architecture.md`; **keep immune and stromal panels separate** (e.g. `results/adata.immune.normalized.p3.h5ad`, `results/adata.stromal.normalized.p3.h5ad`, or equivalent).

---

## 3. Phase Alignment (sc_tools Phasing Scheme)

| Phase | lymph_dlbcl Status | Notes |
|-------|--------------------|--------|
| **1** | Pending | Data ingestion: CSVs (and any mcd/txt) → AnnData; QC. |
| **2** | Pending | Metadata attachment (sample, clinical). |
| **3** | Pending | Preprocessing, normalization, clustering (no cell typing in Phase 3). |
| **3.5b** | Pending | Gene scoring, automated cell typing, optional deconvolution. |
| **4** | Pending | Manual cell typing (if needed). |
| **5–7** | Pending | Downstream biology; optional meta analysis. |

Entry point for this project: start at **planning and inventory** (below), then Phase 1 once CSVs and descriptors are mapped.

---

## 4. Planning Roadmap (Current Priority)

**Priority:** Reuse the previous analyst’s **already-processed data** when the inventory shows it exists and is suitable. **Process raw data only if** we cannot find or use that processed output.

### Step 1: Remote file inventory and directory-structure dataframe (first deliverable)

1. **Obtain remote file listing (no bulk download):**
   - On remote (e.g. `ssh cayuga`): run `find` under `/home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2` for all `*.csv` (and optionally other processed outputs: e.g. `*.h5ad`, `*.rds`, `*.h5`); capture relative path, size, and optionally directory depth.
   - Save raw listing (e.g. `metadata/remote_csv_listing.txt` or equivalent) if useful for reproducibility.

2. **Build a dataframe of CSV (and other key) files with descriptors:**
   - Columns: at least `relative_path`, `remote_full_path`, `size_bytes`, `parent_dir`, `filename`; add **descriptor** and **inferred_role** by cross-referencing:
     - **Project tree:** folder names and path structure (e.g. immune vs stromal preprocessing, stroma_spatial, total_tumor, clinical_analysis).
     - **Notebooks:** scan `notebooks/dlbcl_notebooks/DLBCLv2/**/*.ipynb` for path references to infer which file is used where (e.g. immune-panel merged counts, stromal-panel cell IDs, clinical metadata, **or already-processed objects**).
   - Human-in-loop: review and fill descriptor and `inferred_role`; **flag candidates for “already processed”** (e.g. merged expression matrices, cluster annotations, Seurat/AnnData intermediates).

3. **Store the inventory in project metadata (two artifacts for traceability):**
   - **`metadata/remote_file_inventory.csv`** — Script output (Step 2a); reproducible; `detailed_descriptor` empty.
   - **`metadata/remote_file_inventory_annotated.csv`** — Same columns with **detailed_descriptor** filled (Step 2b: human or AI, or via `scripts/annotate_inventory_descriptors.py`). Use the annotated file for deciding what to download and for documentation.

4. **Document in Journal.md:** decisions on column names, descriptor taxonomy, and any script used to generate the inventory.

### Step 2 (priority): Identify required files and reuse existing processed data

- **What we need to download:** Expression matrices with **cells as observations** (rows). Priority: **(1) Normalized** expression matrix per panel — we will run the sc_tools pipeline on normalized data. **(2) Raw** expression matrix per panel — if we can reassemble or identify raw cell-by-marker tables on the remote, we should download those too so we have both raw and normalized available (e.g. for QC, backup `adata.raw`, or re-normalization). From the inventory and notebook review, identify which remote files correspond to these.
- **Identify required files:** Mark in the inventory (or `metadata/files_to_download.csv`) which paths are: (a) **normalized** expression (cells = observations) for immune and stromal panels; (b) **raw** expression (cells = observations) for each panel when available. Include any cell metadata or sample IDs needed to build AnnData.
- **Download required files:** Download only those identified files to the local project (e.g. `data/downloaded/`). Use a documented method (e.g. `rsync`, `scp`) and keep `metadata/download_manifest.csv` for traceability.
- **Build AnnData:** Use normalized data as the primary input for the sc_tools pipeline. Where raw expression was downloaded, reassemble into AnnData and store raw counts (e.g. in `adata.X` or `adata.layers['raw']`) so we retain both; otherwise build from normalized only.
- **If suitable processed data is found:** Map to our checkpoint names; write `results/adata.immune.*.h5ad` and `results/adata.stromal.*.h5ad` (and shared metadata as needed). Skip re-processing from raw for that panel when normalized (and optionally raw) are already available.
- **If not found (or only partially):** Proceed to Step 3 for the missing panel(s) or samples.

### Step 3 (fallback): Process raw data only when needed

- **Only if** we cannot find or use existing processed data: scripted ingestion from raw CSVs (and any mcd/txt) into AnnData; QC; write `results/adata.immune.raw.p1.h5ad` and/or `results/adata.stromal.raw.p1.h5ad` and then run preprocessing. Keep immune and stromal pipelines separate.

---

## 5. Completed Tasks

- [x] Create project layout: `data/`, `figures/`, `metadata/`, `scripts/`, `results/`, `outputs/`, and planning docs.
- [x] Step 1: Remote file listing obtained (Option B); inventory script run; **annotated inventory** produced and stored as **`metadata/remote_file_inventory_annotated.csv`** (script output kept as `metadata/remote_file_inventory.csv` for traceability).
- [ ] Step 2: Identify required files from annotated inventory; download; build AnnData for immune and stromal panels.

---

## 6. Active Tasks (Roadmap)

- [x] **Step 1.1:** Remote listing obtained (Option B); stored in `metadata/remote_csv_listing.txt`.
- [x] **Step 1.2:** Inventory script run; notebook references and heuristics in `metadata/remote_file_inventory.csv`.
- [x] **Step 1.3:** **Detailed descriptors** filled via `scripts/annotate_inventory_descriptors.py`; output stored as **`metadata/remote_file_inventory_annotated.csv`** (original script output unchanged for traceability).
- [ ] **Step 2 (priority):** Using **`metadata/remote_file_inventory_annotated.csv`**, identify **required** files: (1) normalized expression (cells as observations) per panel; (2) raw expression per panel if available; (3) cell metadata/sample IDs. Record in `metadata/files_to_download.csv`; download to `data/downloaded/`; keep `metadata/download_manifest.csv`; build AnnData into `results/adata.immune.*.h5ad` and `results/adata.stromal.*.h5ad`.
- [ ] **Step 3 (fallback):** For panel(s) with no suitable processed data, run Phase 1 from raw; keep immune and stromal separate.

---

## 7. Traceability and Reproducibility

- **Scripts:** All steps that can be automated (inventory build, notebook scan, download manifest) live under `scripts/` with clear inputs and outputs. Each script logs or writes to `metadata/` or `outputs/` so the pipeline is reproducible.
- **Documented commands:** Remote `find` and download commands (e.g. `rsync`/`scp`) are recorded in this repo (e.g. in `scripts/README.md` or a runbook in `metadata/`) so anyone can re-run the listing and download steps.
- **Manifests:** Keep a manifest of downloaded files (remote path, local path, date, checksum optional) in `metadata/` (e.g. `metadata/download_manifest.csv`) after downloading required files.
- **Inventory:** Script output = `metadata/remote_file_inventory.csv` (reproducible). Annotated version with natural-language descriptors = **`metadata/remote_file_inventory_annotated.csv`** (separate file for traceability; use for downstream decisions).
- **Journal:** Log material decisions (which files were chosen as “required”, download date, script versions) in `Journal.md`.

---

## 8. Blockers and Sanity Checks

- **Blocker:** Do not assume all data is local; plan around remote path access (SSH + find/rsync/scp as needed).
- **Sanity:** Descriptor and inferred_role should be consistent and reviewable so we know which files feed which panel (immune vs stromal) and whether to reuse as processed or ingest as raw.

---

## 9. Technical Decision Log (Reference this project's Journal.md)

- **Inventory columns:** `relative_path`, `remote_full_path`, `size_bytes`, `parent_dir`, `filename`, `descriptor`, `inferred_role` (to be refined in Step 1 and logged in Journal.md).
- **Descriptor source:** Project tree + notebook path references; human review for final labels. **Panels:** immune and stromal; keep separate. **Priority:** reuse existing processed data; process raw only if needed.
- **Required download content:** Normalized expression matrix (cells = observations) per panel — primary; run sc_tools pipeline on normalized. Raw expression matrix (cells = observations) per panel — download when we can reassemble/identify it so we have both raw and normalized (e.g. for adata.raw or re-normalization).
