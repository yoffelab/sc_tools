# Research Journal & Decision Log: lymph_dlbcl (Two-Panel IMC DLBCL)

This journal documents **project-specific** technical and analytical decisions. Repository-level decisions are in the root `Journal.md`.

**Project:** `projects/imc/lymph_dlbcl`  
**Reference:** Root `Mission.md` (toolkit); this project's `Mission.md` (study aims and roadmap).

---

## Log Entries

### [2025-02-12] – Project setup and planning phase

- **Action:** Created project layout and planning documents for the two-panel IMC DLBCL project (lymph_dlbcl). Added standard directories (`data/`, `figures/`, `metadata/`, `scripts/`, `results/`, `outputs/`) and three core docs: `Mission.md`, `Architecture.md`, `Journal.md`.
- **Context:** Picking up work from a prior analyst; data lives remotely at `/home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2`. Processing notebooks are already present under `notebooks/dlbcl_notebooks/DLBCLv2/` (immune and stromal panels, preprocessing, merged analyses, clinical, spatial).
- **Decisions:**
  - **Remote path:** Use backup path above as canonical source; investigate without downloading all data first.
  - **First deliverable:** Build a dataframe of all CSV (and processed) files under the remote path with descriptors; store in `metadata/remote_csv_inventory.csv`; flag processed vs raw.
  - **Inventory columns (planned):** `relative_path`, `remote_full_path`, `size_bytes`, `parent_dir`, `filename`, `descriptor`, `inferred_role`. Refinements to be logged here as Step 1 proceeds.
- **Next:** Run remote `find`; parse notebooks; produce inventory; then prioritize identifying and reusing existing processed data.

### [2025-02-12] – Priority: reuse processed data; panel names (immune, stromal)

- **Action:** Updated Mission.md and Architecture.md to reflect two clarifications from project lead.
- **Decisions:**
  - **Priority:** Reuse the previous analyst’s **already-processed data** when the file inventory shows it exists and is suitable. **Process raw data only if** we cannot find or use that processed output. Step 2 is now “identify and reuse existing processed data”; Step 3 is “process raw only when needed.”
  - **Panel names:** The two panels are **immune** and **stromal** (not stroma/T cell). Keep them **separate** in all checkpoints (e.g. `results/adata.immune.normalized.p3.h5ad`, `results/adata.stromal.normalized.p3.h5ad`).
- **Rationale:** Avoid re-processing when prior work can be picked up; consistent panel naming and separation for downstream analysis.

### [2025-02-12] – Identify required files, download, and traceability

- **Action:** Updated Mission.md and Architecture.md; added reproducible workflow and scripts.
- **Decisions:**
  - **Identify and download:** After building the inventory, we identify which files are **required** (processed, normalized, or filtered) to build AnnData; record them in `metadata/files_to_download.csv`; **download** only those to `data/downloaded/`; keep `metadata/download_manifest.csv` (remote_path, local_path, date) for traceability. Then build AnnData from the downloaded files.
  - **Traceability and reproducibility:** (1) All automatable steps live in `scripts/` with clear inputs and outputs. (2) Remote `find` and download commands are documented in `scripts/README.md` so the listing and download steps can be re-run. (3) Manifests (listing, inventory, files_to_download, download_manifest) are stored in `metadata/`. (4) Material decisions are logged in Journal.md.
- **Artifacts:** Added `scripts/build_remote_inventory.py` (reads `metadata/remote_csv_listing.txt`, optionally scans notebooks for path references, writes `metadata/remote_csv_inventory.csv`); `scripts/README.md` (runbook with exact remote commands, download steps, manifest recording); `metadata/remote_csv_listing_sample.txt` for testing; `data/downloaded/` for required files. Verified inventory script runs with sample listing (with and without `--skip-notebooks`).

### [2025-02-12] – Required files: normalized and raw expression (cells as observations)

- **Action:** Specified in Mission.md what “required files” means for download and pipeline use.
- **Decisions:**
  - **Primary:** Download **normalized** expression matrix with **cells as observations** (rows) per panel (immune, stromal). The sc_tools pipeline will run on this normalized data.
  - **Secondary:** Where we can **reassemble or identify raw** expression matrix (cells as observations) on the remote, download that too so we have both raw and normalized. When building AnnData, store raw counts (e.g. in `adata.X` or `adata.layers['raw']`) when available so we retain both; supports QC, `adata.raw`, or re-normalization.
- **Rationale:** Clear definition of “required” for files_to_download; pipeline is normalized-first with optional raw preserved when possible.

### [2025-02-12] – Annotated file inventory from notebooks (remote_file_inventory.csv)

- **Action:** Refined `scripts/build_remote_inventory.py` to scan all ipynb files, extract context snippets where each listed file is referenced, and write an annotated dataframe to `metadata/remote_file_inventory.csv`.
- **Columns:** relative_path, remote_full_path, size_bytes, parent_dir, filename, descriptor, **data_type** (raw | normalized | processed | metadata | unknown), **analysis_stage** (processing | downstream | clinical | visualization | metadata | other), **detailed_annotation** (notebook cell snippets), **notebook_refs**, **needs_download** (yes | no | review). Listing can include .csv, .h5ad, .rds, .h5 (Option B from runbook).
- **Logic:** Paths in notebooks are normalized and matched to the listing by full path or basename. data_type and analysis_stage are inferred from path/keywords and notebook names; needs_download prioritizes yes for normalized/processed expression and metadata (cell IDs, clinical), then no for clear outputs (write.csv, images), else review.
- **Output:** Script writes to `metadata/remote_file_inventory.csv`. README and Mission/Architecture updated to reference this file and column set.

### [2025-02-12] – Two-step inventory: extract (script) vs annotate (human or AI)

- **Action:** Split the inventory procedure into Step 2a (reproducible Python script) and Step 2b (natural-language annotation of what each file is).
- **Decisions:**
  - **Step 2a (script):** Extracts features from the listing and notebooks (paths, snippets, data_type, analysis_stage, needs_download) and writes `metadata/remote_file_inventory.csv`. A new column **detailed_descriptor** is left empty by the script.
  - **Step 2b (human or AI):** The annotator (human or AI agent) uses `detailed_annotation` and `notebook_refs` to fill **detailed_descriptor** with a short, natural-language summary of what each file contains and how it is used. This is better done by an agent that can read and summarize notebook context than by a static script.
  - **Reproducibility:** Re-running the script overwrites the CSV; use `--merge-descriptors metadata/remote_file_inventory.csv` to copy existing detailed_descriptor values for matching paths so annotations are preserved.
- **Rationale:** Keeps extraction reproducible while allowing flexible, understandable descriptions that require interpreting notebook context.

### [2025-02-12] – Annotated inventory in separate file for traceability

- **Action:** Filled the detailed_descriptor column for all 1,668 rows and stored the result in a **separate file** so the original script output remains unchanged.
- **Decisions:**
  - **Script:** Added `scripts/annotate_inventory_descriptors.py`; it reads `metadata/remote_file_inventory.csv` and writes **`metadata/remote_file_inventory_annotated.csv`** with `detailed_descriptor` filled using rule-based summaries from descriptor, data_type, analysis_stage, detailed_annotation, and path/filename.
  - **Traceability:** `remote_file_inventory.csv` = reproducible script output (unchanged). `remote_file_inventory_annotated.csv` = annotated version for downstream use (what to download, documentation).
- **Roadmap (Mission.md):** Step 1.1–1.3 marked complete; current priority is Step 2 (identify required files from annotated inventory, download, build AnnData).

### [2025-02-12] – RDS notebook usage metadata (rds_notebook_usage.csv, notebook_rds_summary.csv)

- **Action:** Added `scripts/analyze_rds_notebook_usage.py` to scan all notebooks for readRDS/saveRDS and path references to the 47 downloaded RDS files; produced two metadata summaries.
- **Outputs:** (1) **`metadata/rds_notebook_usage.csv`** — per RDS: inferred_content (e.g. Stromal panel S1: Seurat object B-cell subset), read_by_notebooks, written_by_notebooks, read/write snippets, summary. (2) **`metadata/notebook_rds_summary.csv`** — per notebook: rds_read, rds_written, rds_count. Enables traceability from RDS to notebooks and from notebooks to RDS, similar to the CSV inventory.

### [2026-03-04] – Manuscript reproduction pipeline: full Snakemake + 25 scripts

- **Action:** Implemented complete reproducible pipeline for DLBCL IMC manuscript (v7.4) figure reproduction. Hybrid approach: reuse 47 existing Seurat-converted h5ad objects for cell typing/labels; consolidate into sc_tools checkpoint format.
- **Scripts created (25 total):**
  - Phase 0: `download_clinical_metadata.sh`, `validate_h5ad_objects.py`
  - Phase 1: `build_panel_adata.py`, `attach_clinical_metadata.py`, `build_spatial_adata.py`
  - Phase 2-3: `validate_celltypes.py`, `build_lme_classes.py`
  - Phase 4 (13 figures): `fig1_single_cell_atlas.py`, `fig2_lme_classes.py`, `fig3_clinical.py`, `fig4_spatial.py`, `fig5_ml_framework.py`, `supp_fig1_qc_panels.py` through `supp_fig8_extended_survival.py`
  - Phase 6: `validate_figures.py`
- **Snakefile:** Full DAG from download -> validate -> build -> figures -> verify. Rules for each phase; `snakemake --cores 8 all` runs everything.
- **config.yaml:** Project parameters, input paths, LME class definitions, ML params (seed=42, 5-fold CV).
- **Mission.md:** Overhauled with 13-figure plan, implementation phases 0-6, key files reference.
- **Key decisions:**
  - Panels: T2 (immune) + S2 (stromal) as primary (most complete with DLC_code and labels)
  - LME class mapping from k=10 k-means (from DLBCL_case_clustering.ipynb): Cold={0,7,9}, CD206={1,8}, Cytotoxic={2,10}, Stroma={3,4}, T_cell_Regulated={5,6}
  - Runtime: cayuga; data already present at `/home/fs01/juk4007/elementolab/sc_tools/projects/imc/lymph_dlbcl`
  - Stats: BH FDR, lifelines for survival, sklearn RF for ML (Fig 5)
- **Next:** Run `download_clinical_metadata.sh` on cayuga, then `snakemake --cores 8 phase0 phase1`
