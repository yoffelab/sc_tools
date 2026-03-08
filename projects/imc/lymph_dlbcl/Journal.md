# Research Journal & Decision Log: lymph_dlbcl (Two-Panel IMC DLBCL)

This journal documents **project-specific** technical and analytical decisions. Repository-level decisions are in the root `Journal.md`.

**Project:** `projects/imc/lymph_dlbcl`  
**Reference:** Root `Mission.md` (toolkit); this project's `Mission.md` (study aims and roadmap).

---

## Log Entries

### [2026-03-07] — Critical fix: X=zeros in p4 checkpoint; layers['raw'] fallback

- **Root cause:** `adata.immune.celltyped.p4.h5ad` has X matrix = all zeros. Expression data is preserved in `layers['raw']` (Seurat-normalized, centered/scaled values). This caused all heatmaps and marker analyses to show blank/flat output.
- **Discovery path:** Fig 1b heatmap showed flat colors through 5+ iterations. Diagnostic logging revealed marker std = 0.0000 for ALL 40 markers. SSH inspection confirmed `X[0:5,0:5] = all zeros`, `X.min() = X.max() = 0.0`, but `layers['raw']` has real data with negative values (Seurat scale.data).
- **Fix applied to 8 scripts:** Added `data_matrix = adata.layers.get("raw", adata.X) if adata.layers else adata.X` pattern to:
  - `fig1_single_cell_atlas.py` — UMAP (copies layers['raw'] to X for PCA), heatmap, B cell markers
  - `supp_fig1_qc_panels.py` — marker intensities and distributions
  - `supp_fig2_bcell.py` — B cell heatmap
  - `supp_fig3_tcell_myeloid.py` — T cell/myeloid heatmap
  - `supp_fig4_vessel.py` — endothelial violin plots
  - `build_lme_classes.py` — TME feature profiles for cluster->LME mapping
  - `validate_celltypes.py` — validation heatmap
  - `validate_h5ad_objects.py` — data quality metrics (reports x_source)
- **Additional fig1 fixes (from prior session, finalized here):**
  - `figure_config.py`: Added `MORPHOLOGICAL_PATTERNS`, `filter_protein_markers()`, `normalize_celltype()`, `build_celltype_palette(normalize=False)` for distinct colors per subtype, `is_bcell_label()`
  - Fig 1a UMAP: 19 distinct colors via tab20+tab20b fallback
  - Fig 1b heatmap: Only T*/M* prefixed subtypes (excludes bcell/stroma/other), z-scored per marker, 40 markers x 16 cell types — shows clear biological patterns (FoxP3 in Tregs, CD68 in macrophages, GrB in cytotoxic)
  - Fig 1c prevalence: Immune panel only, excludes bcell/Unknown/other
- **Fig 2-4 safe:** These scripts only use `adata.obs` metadata, not X matrix.
- **Verified:** Fig1 UMAP, heatmap, prevalence all produce correct output. Supp fig 1 distributions show proper marker intensities.

### [2026-03-07] — Plan B: IMC reprocessing infrastructure (Steps 1-4 of 5)

- **Action:** Implemented Plan B pipeline infrastructure for raw IMC reprocessing from Hyperion MCD files on cayuga. Steps 1-4 complete; Step 5 (SLURM submission) blocked on cayuga VPN access.
- **Files created/modified:**
  - `scripts/generate_panel_csv.py` — extracts protein names from h5ad var_names (`Protein-mem`/`-nuc` pattern), filters morphological features, writes imctools-compatible panel CSVs with empty Metal_Tag (to fill from MCD headers). Immune: 38 channels; stromal: 38 channels.
  - `scripts/generate_phase0_manifests.py` — reads DLC codes from immune h5ad (`DLC_code` column; 84 unique), falls back to clinical TSV for stromal panel (349 cases; stromal h5ad has no patient-level tracking in `orig.ident`). Writes `metadata/phase0/batch1_immune.tsv`, `batch1_stromal.tsv`, `all_samples.tsv`. mcd_file paths are placeholders (`{backup_dir}/{sample_id}.mcd`).
  - `metadata/panel_immune_t2.csv`, `metadata/panel_stromal_s2.csv` — generated (Metal_Tag empty).
  - `metadata/phase0/batch1_immune.tsv` (84 rows), `batch1_stromal.tsv` (349 rows), `all_samples.tsv` (433 rows) — generated with placeholder paths.
  - `config.yaml` — added `phase0a` block (enabled: false, raw_data_dir, output_dir, panel/manifest paths, pipeline_script).
  - `Snakefile` — added `phase0a_immune`, `phase0a_stromal`, `run_imc_pipeline` rules gated by `_PHASE0A_ENABLED` flag; reads manifests only when enabled.
- **Key findings from data exploration:**
  - Immune h5ad (`t2_SO_seurat.h5ad`): 138 unique DLC codes in `obs['DLC_code']` (84 match `DLC\d+` format, rest are CTMA121/0); 49 var_names in `Protein-compartment` format (CellProfiler output, NOT `Protein(IsotopeTag)` format).
  - Stromal h5ad (`S1_seurat_SO.h5ad`): no DLC codes in any obs column; `orig.ident` = integer cluster IDs; `full_index` = cell row indices. Stromal patient tracking requires `S2_seurat_cellid.csv` from cayuga (not yet downloaded) or path sheet CSVs.
  - var_names do NOT contain isotope tags — Metal_Tag column must be filled manually from MCD file channel headers or original acquisition panel sheet.
- **cayuga status:** Not reachable from macOS without Cornell VPN (connection timeout at `cayuga-login1.cac.cornell.edu`). brb is reachable but does not mount cayuga filesystems.
- **Blocking items for Step 5:**
  1. Connect to cayuga (VPN) to run discovery commands from Plan Step 1.
  2. Verify actual MCD file paths (may differ from `{backup_dir}/{sample_id}.mcd` pattern; `3.28.23.DLBCL_image_sheet.csv` on cayuga has authoritative paths).
  3. Fill `Metal_Tag` in panel CSVs from MCD headers or `6.22.22.file_path_sheet_s1/s2.csv`.
  4. Verify imctools/CellProfiler availability on cayuga (`conda activate sc_tools && python -c "import imctools"`).
  5. Confirm `~/elementolab/imc/run_pipeline.py` path for `pipeline_script`.
  6. Run Snakemake dry-run and submit.
- **Decision:** Stromal manifest uses clinical TSV fallback (--use-clinical-for-stromal flag) since stromal h5ad has no patient tracking. Once path sheets are downloaded from cayuga, re-run manifest generation with `--path-sheet-stromal` to populate real mcd_file paths and sample IDs.

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

### [2026-03-06] - Pipeline v2: Data fixes + figure quality overhaul

- **Action:** Major overhaul addressing missing data, incorrect figures, and codifying figure production principles.
- **Problems identified from run8 (46 PDFs):**
  - Fig 1 missing (no PDFs generated)
  - Stromal panel had only 1 sample (orig.ident not mapped to DLC_code)
  - No survival curves (column name mismatch)
  - LME proportions wrong (missing Cold class; 4 classes instead of 5)
  - Heatmaps not z-scored, inconsistent colors, no statistical annotations
- **Changes:**
  - **skills.md Section 12.8:** Added "Figure Intent, Insight, and Readability" — separates general figure principles (every figure needs a claim, direction-of-effect validation) from reproduction-specific rules (align with manuscript text, quantitative validation). Key principle: if a bar plot shows a statistical comparison, the direction must match the stated insight.
  - **figure_config.py:** Shared infrastructure — LME_COLORS (5 classes), LME_ORDER (manuscript order), CELLTYPE_COLORS (12 types, Okabe-Ito), COO_COLORS (GCB/ABC), FIGURE_RC_PARAMS (Arial 7pt, 300 DPI, PDF type 42). All figure scripts import from here.
  - **attach_clinical_metadata.py:** Rewrote to use metadata/DLC380_clinical.tsv directly (348 cases, tab-separated). No more download dependency. Standardized column names: OS_time, OS_event, PFS_time, PFS_event, COO, age, sex, IPI, stage. Filters FINAL_COHORT=YES (332 cases).
  - **build_lme_classes.py:** Three-tier LME assignment: (1) patient_tme CSV with direct names, (2) TME cluster CSV with k=10 mapping, (3) k-means on abundance fallback. Validates proportions against manuscript (Cold 35.1%, Stromal 21.3%, Cytotoxic 20.7%, T cell Regulated 14.6%, CD206 Enriched 8.2%). Warns if any class off by >5%.
  - **build_panel_adata.py:** Added barcode parsing fallback for stromal panel — regex extracts DLC_code from cell barcodes when orig.ident gives only 1 sample.
  - **Fig 1:** Added panels c (non-B cell prevalence), d (B cell markers), e (COO distribution). Z-scored marker heatmap. CELLTYPE_COLORS palette.
  - **Fig 2:** Heatmap z-scored per cell type across LME classes (not raw). LME_ORDER for consistent row ordering. BH-corrected 1-vs-rest Wilcoxon on violins.
  - **Fig 3:** Reads DLC380_clinical.tsv directly. Multivariate log-rank test. Cox forest with HR/CI/p-values (CD206 Enriched as reference). COO_COLORS palette.
  - **Fig 4:** Community heatmap z-scored per community. Shannon entropy bar chart for diversity. Community x LME enrichment bubble plot (orange=enriched, purple=depleted). Gradient boosted model (not RF) for prediction.
  - **Fig 5:** Moved ML analysis to `figures/manuscript/extended/ml_classification/` (manuscript Fig 5 is CosMx, deferred until data available).
  - **All supp figs:** Updated to import figure_config and apply_figure_style(). Removed duplicate LME_COLORS definitions.
  - All scripts pass ruff lint.
- **Key decisions:**
  - DLC380_clinical.tsv is the canonical clinical metadata source (already in metadata/, no download needed)
  - Fig 5 is CosMx (deferred) — ML analysis is extended data
  - Gradient boosted model for LME prediction (matching manuscript, not RF)
  - Direction-of-effect validation is now a codified principle in skills.md
- **Next:** Upload to cayuga and run: Phase 0 (download metadata CSVs), Phase 1 (build panel adata), Phase 3 (LME assignment), Phase 4 (figures). Then rsync figures back for visual QA.
