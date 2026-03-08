# Mission: sc_tools Toolkit and Pipeline

**Scope:** Repository-level and generalizable goals. Project-specific objectives (e.g. TLS, macrophage analysis) live in each project's `Mission.md` under `projects/<data_type>/<project_name>/Mission.md`.

**Last Updated:** 2026-03-06

---

## 1. Toolkit and Pipeline (General)

- **Pipeline phases (0–7):** Non-linear workflow with human-in-loop steps. See `README.md` for the diagram and phase summary.
- **sc_tools:** Reusable package (`sc_tools.pl`, `sc_tools.tl`, `sc_tools.qc`, etc.) for QC, plotting, testing, colocalization, I/O. Keep generic; project-specific logic stays in project scripts.
- **Reproducibility:** Makefile is project-aware (`PROJECT ?= projects/visium/ggo_visium`). Each project has `data/`, `figures/`, `metadata/`, `scripts/`, `results/`, `outputs/` under `projects/<platform>/<project_name>/`.
- **Standards:** All projects follow `skills.md` (statistics, significance bars, FDR). Documentation avoids apostrophes.

### Phasing scheme (sc_tools)

| Phase | Name | Checkpoint | Required Data | QC Report |
|-------|------|------------|---------------|-----------|
| **ingest_raw** (0a) | Platform tools (HPC) | `data/{sample_id}/outs/` | Platform-specific raw output | |
| **ingest_load** (0b) | Load per-sample AnnData | `data/{sample_id}/adata.ingested.h5ad` | `obs[sample, library_id, raw_data_dir]`, `obsm[spatial]`, `X` raw counts | |
| **qc_filter** (1) | QC and Concatenation | `results/adata.filtered.h5ad` | `obs[sample, raw_data_dir]`, `obsm[spatial]`, `X` raw, concatenated | `pre_filter_qc.html` |
| **metadata_attach** (2) | Metadata Attachment | `results/adata.annotated.h5ad` | + clinical columns in `obs` | `post_filter_qc.html` |
| **preprocess** (3) | Preprocessing | `results/adata.normalized.h5ad` | `obsm[X_scvi or embedding]`, `obs[leiden]`, `adata.raw` backed up | `post_integration_qc.html` |
| **demographics** (3.5) | Demographics | Figure 1 | Cohort metadata from preprocess | |
| **scoring** (3.5b) | Gene scoring / Auto cell typing | `results/adata.scored.h5ad` | `obsm[signature_score, signature_score_z]`, `uns[signature_score_report]` | |
| **celltype_manual** (4) | Manual Cell Typing | `results/adata.celltyped.h5ad` | + `obs[celltype, celltype_broad]` | `post_celltyping_qc.html` |
| **biology** (5) | Downstream Biology | `figures/manuscript/` | Reads from scoring or celltype_manual checkpoint | |
| **meta_analysis** (6-7) | Meta Analysis | `results/adata.{level}.{feature}.h5ad` | `obs` indexed by roi/patient; `X` = aggregated feature | |

All QC reports are date-versioned (`YYYYMMDD`) and saved to `figures/QC/`. Full validation contracts: Architecture.md Section 2.2.

**Entry points:** Start at Phase 3 (preprocessed AnnData), Phase 3.5b (clustered AnnData), or Phase 4 (phenotyped AnnData). See README diagram for START HERE conditions.

**Checkpoint nomenclature:** Standard filenames are defined in **Architecture.md** (Section 2.1): `adata.ingested.h5ad` (per-sample), `adata.filtered.h5ad` (concatenated QC), `adata.annotated.h5ad`, `adata.normalized.h5ad`, `adata.scored.h5ad`, `adata.celltyped.h5ad`, and `adata.{level}.{feature}.h5ad` for aggregated. Legacy p-code names (`adata.raw.p1.h5ad`, `adata.p0.h5ad`, `adata.raw.h5ad`, etc.) are accepted during transition but deprecated.

---

## 2. General Analysis — To-Do and Implementation Status

All paths below are project-specific: `projects/<platform>/<project_name>/...`. Nothing lives at repo root for `metadata/`, `results/`, or `figures/`.

### Phase 0: Upstream Raw Data Processing

#### Phase 0a: Platform tools (HPC)

- [x] **Batch manifest system:** `sc_tools.ingest.config` — `load_batch_manifest()`, `collect_all_batches()`, `validate_manifest()`. Per-batch TSVs under `metadata/phase0/`; auto-collects into `all_samples.tsv`.
- [x] **Space Ranger command builder:** `sc_tools.ingest.spaceranger` — `build_spaceranger_count_cmd()`, `build_batch_commands()`. Supports Visium (--image) and Visium HD (--cytaimage).
- [x] **Xenium Ranger command builder:** `sc_tools.ingest.xenium` — `build_xenium_ranger_cmd()`.
- [x] **IMC pipeline command builder:** `sc_tools.ingest.imc` — `build_imc_pipeline_cmd()`.
- [x] **Snakemake Phase 0a rules:** `spaceranger_count` and `phase0` target in robin, ggo_visium, and create_project.sh template. Driven by `metadata/phase0/all_samples.tsv`.

#### Phase 0b: Load into AnnData / SpatialData

- [x] **Modality loaders (implemented):** `sc_tools.ingest.loaders` — `load_visium_sample()`, `load_visium_hd_sample()`, `load_visium_hd_cell_sample()`, `load_xenium_sample()`, `load_imc_sample()`. Each sets `obs['sample']`, `obs['library_id']`, `obs['raw_data_dir']`, `obsm['spatial']`. Output: `data/{sample_id}/adata.ingested.h5ad`.
- [ ] **CosMx loader (deprioritized):** `load_cosmx_sample()` — reads flat CSV/Parquet files or RDS (via rpy2+anndata2ri) from NanoString/AtoMx output. Sets spatial coordinates from cell centroid (x, y in microns). Batch TSV schema: `sample_id`, `cosmx_dir`, `batch`. Three panel tiers with different output characteristics:
  - **CosMx 1k:** ~1,000-plex targeted panel; flat CSV/Parquet export from AtoMx.
  - **CosMx 6k:** ~6,000-plex panel; larger expression matrices; same flat file format.
  - **CosMx full_library:** Whole-transcriptome (~18k genes); significantly larger files; may require chunked loading or backed AnnData.
- [x] **visium_hd_cell modality:** SpaceRanger 4 cell segmentation support. `load_visium_hd_cell_sample()` reads from `outs/cell_segmentation/`. Cell-level QC thresholds (like xenium). Routes to xenium recipe in `preprocess()`. `create_project.sh` accepts `visium_hd_cell` as valid data type.
- [x] **Concatenation:** `concat_samples()` merges per-sample AnnData with `calculate_qc_metrics()` applied.
- [ ] **IMC end-to-end ingestion pipeline (priority):** Wire `build_imc_pipeline_cmd()` + `load_imc_sample()` into a complete Snakemake workflow for IMC projects. Includes: Phase 0a Snakemake rules (run ElementoLab pipeline per MCD on HPC → `processed/{sample}/`), Phase 0b rules (load each sample → `data/{sample_id}/adata.ingested.h5ad`), and batch manifest templates for IMC (`metadata/phase0/`). Must be done **before** the generic checkpoint script below. Target projects: lymph_dlbcl, ggo-imc.
- [x] **IMC image loading (Phase 0b):** `load_imc_sample(load_images=True, panel_csv=...)` reads the per-ROI `{roi_id}_full.tiff` multi-channel stack (C, H, W) and companion `{roi_id}_full.csv` channel index file from `processed/{sample}/tiffs/`. Builds arcsinh-normalized full stack + RGB composite (PanCK=R, CD3=G, DNA1=B defaults). Stores in `adata.uns['spatial'][sample_id]` (squidpy/scanpy-compatible). `IMCPanelMapper` resolves names from `MarkerName(IsotopeTag)` format (exact, case-insensitive, isotope tag, partial); reads both `*_full.csv` (per-ROI channel index) and `channel_labels.csv` (project panel with full/ilastik flags). Optional mask loading (`load_mask=True`) for `*_full_mask.tiff`. 18 new unit tests; 74 total pass, lint clean.
- [x] **H&E image loading (any modality):** `load_he_image(he_path, library_id, adata)` injects any TIFF (RGB, grayscale, or CHW format) into `adata.uns['spatial'][library_id]['images']`. Handles grayscale-to-RGB, uint16/float-to-uint8, channel-first-to-HWC, and downsample. 5 unit tests.
- [x] **IMC spatial plotting:** `sc_tools.pl.plot_imc_composite()` — same API as `plot_spatial_plain_he`; renders `images['hires']` RGB composite with channel label title. `sc_tools.pl.plot_imc_channel()` — renders single channel from `images['full']` (C, H, W) stack with inferno colormap and percentile-based vmax. Both exported in `sc_tools.pl.__all__`.
- [x] **Phase 0b checkpoint script (`scripts/ingest.py`):** Generic CLI reads batch manifest, dispatches to modality loaders (`_get_loader()`), saves per-sample `data/{sample_id}/adata.ingested.h5ad` (`--save-per-sample`), concatenates, saves final output. 5 tests in `test_ingest_script.py`.
- [ ] **Phase 0b Snakemake rule:** `adata_ingested` rule with sentinel `data/{sample_id}/.adata.ingested.done`. Wire `scripts/ingest.py` into project Snakefiles.
- [ ] **SpatialData (optional):** `data/{sample_id}/spatialdata.zarr` for Visium HD and Xenium when full image pyramids / subcellular coords needed. Loader via `spatialdata-io`.

### Checkpoint Validation

- [x] **`sc_tools.validate`:** `validate_p1` through `validate_p4`, `validate_checkpoint()` dispatcher, `validate_file()` convenience. Auto-fix: `obs['batch']` to `obs['raw_data_dir']` (p1). 30 unit tests.
- [x] **CLI:** `scripts/validate_checkpoint.py` — `--phase`, `--fix`, `--warn-only`. Exit code 1 on failure.
- [x] **Snakemake validation sentinels:** `results/.adata.{name}.validated` sentinel files in robin, ggo_visium, and create_project.sh template. Downstream rules depend on sentinels.

### Phase 1: QC and Concatenation

**Input:** Per-sample `data/{sample_id}/adata.ingested.h5ad` produced in Phase 0b.

- [x] **Load Phase 0 checkpoints:** Read all per-sample `adata.ingested.h5ad` files listed in `metadata/phase0/all_samples.tsv`. Implemented in `scripts/ingest.py`.
- [x] **Per-sample QC:** `sc_tools.qc.filter_spots()` (modality-aware thresholds); `compute_sample_metrics()`; `classify_samples()` (absolute + MAD outlier thresholds).
- [x] **Concatenation:** `concat_samples()` across all passing samples → `results/adata.filtered.h5ad`. Implemented in `sc_tools.ingest`.
- [x] **Required annotations:** `obs['sample']`, `obs['raw_data_dir']`, `obsm['spatial']`; `X` raw counts; no normalization.

#### QC Metrics (sc_tools.qc)
- [x] Implement `calculate_qc_metrics` (wrap scanpy `pp.calculate_qc_metrics`; optional mt/hb patterns).
- [x] Implement `filter_cells`, `filter_genes` (wrap scanpy).
- [x] Implement `highly_variable_genes` (wrap scanpy).
- [x] Implement `spatially_variable_genes` (wrap squidpy `gr.spatial_autocorr`). Sample-specific QC on spots/cells and genes/proteins.
- [x] Two QC versions: pre-normalization and post-normalization (same API; call at raw and after normalize).

#### QC Plots
- [x] 2x2 grid: total_count (total expression), gene/protein expression histogram, log1p version (all platforms). `sc_tools.qc.plots.qc_2x2_grid`; also `sc_tools.pl.qc_2x2_grid`.
- [x] 2x4 pre vs post: left 2x2 pre-filter, right 2x2 post-filter. `sc_tools.qc.plots.qc_2x4_pre_post`.
- [x] Spatial transcriptomics: % mt, % hb per spot (from `calculate_qc_metrics` with mt_pattern/hb_pattern).
- [x] Multipage spatial plot per sample: total_count, log1p, %mt (1x3 subplot) as QC report. `sc_tools.qc.plots.qc_spatial_multipage`; also `sc_tools.pl.qc_spatial_multipage`.
- [x] Violin metrics: n_genes_by_counts, total_counts, pct_counts_mt (optional groupby). `sc_tools.qc.plots.qc_violin_metrics`.
- [x] Scatter: total_counts vs n_genes_by_counts, color pct_counts_mt. `sc_tools.qc.plots.qc_scatter_counts_genes`.
- [x] Highly variable genes plot: mean vs dispersion, HVG highlighted. `sc_tools.qc.plots.plot_highly_variable_genes`.
- [x] Spatially variable genes plot: mean/rank vs Moran's I; per-library when `library_id` present (`spatially_variable_genes_per_library`). Skip SVG only when library_id missing. `sc_tools.qc.plots.plot_spatially_variable_genes`.
- [x] Save under `projects/<platform>/<project>/figures/QC/raw/` (pre) and `figures/QC/post/` (post). Caller passes output_dir or output_path.
- [x] Single QC report script: `scripts/run_qc_report.py` produces all QC figures; Snakemake rule `qc_report` → `figures/QC/qc_report.done`.
- [x] Cross-sample comparison plots: `qc_sample_comparison_bar`, `qc_sample_violin_grouped`, `qc_sample_scatter_matrix` with pass/fail coloring.

#### Sample-Level QC (sc_tools.qc.sample_qc)
- [x] `filter_spots()`: modality-aware spot filtering with conservative defaults (visium/visium_hd/xenium/cosmx/imc).
- [x] `compute_sample_metrics()`: per-sample aggregate QC (n_spots, median counts/genes, %MT stats, n_genes_detected, Space Ranger metrics).
- [x] `classify_samples()`: absolute thresholds + adaptive MAD outlier detection; conservative defaults minimize drops.
- [x] `save_pass_fail_lists()`: write `qc_sample_pass.csv` and `qc_sample_fail.csv`.
- [x] `apply_qc_filter()`: backup + spot filter + sample removal + save.

#### QC HTML Report (sc_tools.qc.report)
- [x] `generate_qc_report()`: self-contained HTML with summary cards, per-sample metrics table, embedded plots. Jinja2 template at `sc_tools/data/qc_report_template.html`. **Deprecated** in favor of three date-versioned reports below.
- [x] Integrated into `scripts/run_qc_report.py` with `--modality`, `--mad-multiplier`, `--thresholds`, `--apply-filter` CLI args.
- [x] 15 new unit tests (28 total in test_qc.py); all pass, lint clean.

#### Date-Versioned QC Reports (restructured — now 4 reports)
- [x] **`generate_pre_filter_report()`**: Phase 1 entry. Output: `figures/QC/pre_filter_qc_YYYYMMDD.html`.
- [x] **`generate_post_filter_report()`**: Phase 1-2 exit. Pre-vs-post comparison, HVG/SVG. Output: `figures/QC/post_filter_qc_YYYYMMDD.html`.
- [x] **`generate_post_integration_report()`**: Phase 3 exit. UMAP grid, cluster distribution, integration benchmark. **Batch score is the primary metric** (celltype labels are preliminary). Output: `figures/QC/post_integration_qc_YYYYMMDD.html`. Accepts `comparison_df` to skip recomputation.
- [x] **`generate_post_celltyping_report()`**: Phase 4 exit. Full bio + batch metrics with validated celltypes. Re-scores `results/tmp/integration_test/{method}.h5ad` if available. Output: `figures/QC/post_celltyping_qc_YYYYMMDD.html`. 3 tests.
- [x] **`generate_all_qc_reports()`**: Orchestrator. `report_utils.py` for shared helpers.
- [x] **`sc_tools.pl.qc_plots`**: `qc_umap_grid()`, `qc_cluster_distribution()`.
- [x] **`bm/integration.py`**: `celltype_key` now optional (None = batch-only metrics). `compare_integrations()` with `comparison_df` pass-through.
- [x] **CLI**: `--report pre_filter|post_filter|post_integration|post_celltyping|all`. New args: `--adata-integrated`, `--batch-key`, `--celltype-key`, `--embedding-keys`, `--segmentation-masks-dir`.
- [x] **Snakemake**: Three QC rules replace one (ggo_visium + create_project.sh). Optional segmentation scoring. Fourth rule for post-celltyping pending.
- [x] 12 new tests; 58 total pass, lint clean.

- [ ] **Future:** MA plots (per gene/protein across samples).

---

### Phase 2: Metadata Attachment (Human-in-Loop)

- [ ] **Bypass file:** `projects/<platform>/<project>/metadata/sample_metadata.csv` or `.xlsx` mapping `sample` → clinical columns. sc_tools guides the join when provided.
- [ ] Without file: human prepares map; cannot skip via automatic pipelining.
- [ ] Must be done as early as possible in the project cycle.
- [ ] Preprocessed projects may skip to Phase 3 or 4 if already annotated.

---

### Phase 3: Preprocessing

- [x] **sc_tools.pp module:** Modality-aware preprocessing with GPU auto-detection (rapids-singlecell fallback to scanpy).
- [x] **Normalization:** `backup_raw`, `normalize_total`, `log_transform`, `scale`, `arcsinh_transform` (IMC), `filter_genes_by_pattern` (MT/RP/HB).
- [x] **Integration:** `run_scvi` (default for Visium), `run_harmony` (PCA-based), `run_cytovi` (IMC protein data). All soft dependencies with install hints.
- [x] **Reduction/Clustering:** `pca`, `neighbors` (auto-detects use_rep), `umap`, `leiden`, `cluster` (convenience wrapper). `run_utag` for spatial-aware clustering (soft dep).
- [x] **Recipes:** `preprocess(adata, modality="visium"|"visium_hd"|"xenium"|"cosmx"|"imc", integration=..., ...)` dispatches to modality-specific pipelines. 38 unit tests, lint clean.
- [x] **Integration methods (sc_tools.pp.integrate):** `run_scvi`, `run_scanvi`, `run_harmony` (direct harmonypy call with PyTorch tensor fix), `run_cytovi`, `run_imc_phenotyping` (log1p → per-ROI z-score → ComBat → BBKNN). ComBat via `scanpy.pp.combat`. All soft deps.
- [x] **Harmony PyTorch fix:** Direct `harmonypy.run_harmony()` call bypassing scanpy wrapper; handles `Z_corr` as torch tensor (`.detach().cpu().numpy()`) or numpy array; validates shape orientation; stores `np.float32` contiguous array.
- [x] **Integration benchmark workflow (`sc_tools.bm.integration`):** `run_full_integration_workflow()` orchestrator: subsample → run all candidate methods → benchmark with batch score as primary metric → select best → apply to full dataset. Stores intermediates in `results/tmp/integration_test/{method}.h5ad`. Records selected method in `results/integration_method.txt`. Helpers: `_stratified_subsample()`, `_resolve_best_method()`, `_apply_integration_method()`. 6 tests.
- [x] **4 QC HTML reports:** Pre-filter (Phase 1), post-filter (Phase 1-2), post-integration (Phase 3, batch-focused), **post-celltyping** (Phase 4, full bio metrics). All 4 implemented. See Architecture.md Section 2.5.
- [x] **`generate_post_celltyping_report()`:** Re-evaluate integration with validated celltypes. Optionally re-score all `results/tmp/integration_test/{method}.h5ad` candidates. Output: `figures/QC/post_celltyping_qc_YYYYMMDD.html`. 3 tests.
- [ ] Post-normalization QC report → `figures/QC/post/`.
- [ ] **No automated cell typing here** (moved to Phase 3.5b).

---

### Phase 3.5: Demographics (Branching)

- [ ] **sc_tools helpers:** piechart, histogram, violinplot, kdeplot, barplot, stacked barplot, scatterplot, correlogram, heatmap.
- [ ] Figure 1 for population-based studies describing the cohort.

---

### Phase 3.5b: Gene Scoring, Automated Cell Typing, Deconvolution

- [x] **Gene signature storage:** Implemented in `sc_tools.tl.score_signature`. Scores stored in **`obsm['signature_score']`** and **`obsm['signature_score_z']`** (single flat key; full-path column names). Report in `uns['signature_score_report']`. Downstream use `get_signature_df(adata)` / `get_signature_columns(adata)` (obsm-first, fallback to obs).
- [x] **Multi-method scoring:** `score_signature(method=...)` supports `"scanpy"` (default), `"ucell"` (rank-based AUC; pyucell soft dep), `"ssgsea"` (gseapy soft dep). Records `uns["scoring_method"]`. 80 unit tests total.
- [x] **Gene set loaders:** `sc_tools.tl.gene_sets` — `load_hallmark()` (50 bundled MSigDB Hallmark sets; stripped from robin metadata; offline), `load_msigdb_json()`, `load_gmt()`, `list_gene_sets()`. Curation: `validate_gene_signatures()`, `merge_gene_signatures()`, `update_gene_symbols()`, `save_gene_signatures()`.
- [x] **Group-level enrichment testing:** `sc_tools.tl.run_ora()` (Fisher exact + BH FDR), `sc_tools.tl.run_gsea_pseudobulk()` (pseudobulk prerank via gseapy; soft dep). Plotting: `sc_tools.pl.plot_gsea_dotplot()`. Optional deps: `pip install sc-tools[geneset]`.
- [x] **Always apply:** Use `merge_gene_signatures(project_sigs, load_hallmark())` before calling `score_signature`. Project JSON + bundled Hallmark combined in one pass.
- [x] **Automated cell typing (cluster → celltype):** `sc_tools.tl.celltype` subpackage. `annotate_celltypes(adata, method=...)` dispatcher with Protocol-based backend registry. Backends: ScType (pure numpy, cluster or per-cell, min_score threshold), custom hierarchical gating (DFS traversal, zero deps), CellTypist (soft dep), ensemble majority/weighted vote. Stubs for scArches, Geneformer, scGPT, SingleR. `apply_celltype_map()` writes `obs[celltype]`, `obs[celltype_broad]`, `uns[celltype_colors]` from JSON map. Bundled marker DBs (`immune_human_rna.json`, `immune_human_protein.json`, `protein_aliases.json`) and gating templates (`immune_broad.yaml`, `lymph_node.yaml`). `pyproject.toml` extras: `[celltyping]`, `[foundation]`. 37 unit tests; 574 total pass.
- [x] **Cell-type deconvolution module:** `sc_tools.tl.deconvolution()` with backend registry (cell2location, tangram, destvi). `extract_reference_profiles()` for memory-optimized Cell2location. Per-library backed loading. Output: `obsm['cell_type_proportions']`, `obs['{method}_argmax']`. Thin wrappers in ggo_visium and robin. 14 unit tests.
- [ ] Phase 3.5b is a separate branch from Phase 3 (parallel to 3.5); connects to Phase 4. Phase 4 is skippable if automated cell typing is adequate.

---

### Phase 4: Manual Cell Typing (Human-in-Loop, Iterative; Skippable)

- [ ] **sc_tools workflow:** Extract `cluster_id` from `adata.obs`; provide JSON template.
- [ ] **JSON format:** `{cluster_id: {celltype_name: "...", total_obs_count: N}}`. Include counts per cluster for user guidance.
- [ ] **cluster_id type:** Match `adata.obs` (string vs int; handle categorical with care to avoid NaN/errors).
- [ ] **Two mappings:** `cluster_id → celltype` (full, e.g. `Epithelial (Ki67/MKI67+)`), `cluster_id → celltype_broad` (without parenthesis, e.g. `Epithelial`).
- [ ] **Per iteration:** matrixplot, UMAP (celltype, celltype_broad, clusters), scatter plots (raw/normalized, gene), signature genes per celltype.
- [ ] Repeat until satisfactory. Save `celltype_map.json` under `projects/<platform>/<project>/metadata/`.
- [x] **Post-celltyping QC report:** `generate_post_celltyping_report()` → `figures/QC/post_celltyping_qc_YYYYMMDD.html`. Purple header distinguishes from other reports. New `post_celltyping_qc_template.html`: celltype composition section (stacked bar + table), Phase 3 vs Phase 4 benchmark table (bio column primary/bold), UMAP grid, celltype distribution per sample, optional marker dotplot. `ranking_rows` bug fixed (was missing). New params: `comparison_df_p3`, `marker_genes`. `qc_celltype_abundance()` added to `sc_tools.pl`. Tab navigation across all 4 reports (inline embed — self-contained HTML). Reports 2 and 3 also embed prior reports as tabs. `_find_latest_report`, `_extract_body_content`, `_extract_head_css`, `_wrap_with_tabs` added to `report_utils.py`. `qc_post_celltyping` Snakemake rule added to robin and ggo_visium. 17 new tests; 574 total pass.

---

### Phase 5: Downstream Biology

- [ ] Uses gene scores (and optionally deconvolution) from Phase 3.5b.
- [ ] Gene set analysis; statistical comparison (1-vs-rest, pairwise per `skills.md`).
- [ ] Spatial/process analysis: colocalization, Moran's I, neighborhood enrichment, ligand-receptor.
- [ ] Visualization: publication-ready figures under project `figures/`.

---

### Phase 6–7: Meta Analysis (Optional)

- [ ] **Phase 6:** Aggregate to ROI level and patient level (mean expression, celltype counts).
- [ ] **Phase 7:** Downstream analysis on aggregated ROI/patient data.

---

## 3. Project-Specific Missions

- **GGO Visium:** `projects/visium/ggo_visium/Mission.md` (aligned with phasing scheme above).
- **Other projects:** Add `Mission.md` and `Journal.md` via `./projects/create_project.sh <project_name> <data_type>`.
- **Entry points:** See Phasing scheme (Section 1). Start at Phase 3 (preprocessed), Phase 3.5b (clustered), or Phase 4 (phenotyped) as appropriate.

---

## 4. Workflow Diagram

See **README.md** (Pipeline Workflow section) for the Mermaid diagram and phase summary.

---

## 5. Testing Strategy

Two-layer testing ensures code compiles, passes tests, and the pipeline runs correctly before and during sc_tools implementation.

| Layer | Location | Purpose |
|-------|----------|---------|
| **Package tests** | `sc_tools/tests/` | Unit tests for sc_tools APIs (pl, tl, qc, etc.) with synthetic fixtures. |
| **Project tests** | `projects/<platform>/<project>/tests/` | Integration/smoke tests for the project pipeline (Makefile, scripts, data flow). |

**Implementation order:**
1. **Project unit tests (ggo_visium)** — Establish baseline; validate Makefile and scripts (Phases 1–5) run with fixtures.
2. **sc_tools unit tests** — Ensure guaranteed behavior of sc_tools functions before use in scripts.
3. **Implement functions** — Refactor scripts to use sc_tools; new code must compile and pass both test layers.

All new code must compile and pass tests. Project scripts that use sc_tools should import from `sc_tools` and live under `projects/<platform>/<project>/scripts/`.

---

## 6. Reorganization & Implementation Status

### Completed

- **Docker + conda + UV:** Dockerfile uses miniconda3 with conda env `sc_tools`; UV installs sc_tools. [project_setup.md](project_setup.md) documents build, run, per-project usage. Robin has run_docker.sh.
- **Journal and Mission-as-todo workflow:** Journal.md at root and per project (lymph_dlbcl, ggo_visium); Mission.md is the todo list; in work mode the agent updates Mission after each prompt. Skill (`.claude/skills/journal-and-mission/`); create_project.sh creates Journal.md for new projects.
- **sc_tools skills as Cursor skill:** Repository root `skills.md` is exposed as Cursor skill `.cursor/skills/sc-tools-skills/SKILL.md`; agent follows skills.md for analysis and coding. Sandbox/local defaults: Docker + Snakemake (documented in skills.md §13 and in the skill).
- **sc_tools package:** `pl/` (spatial, heatmaps, statistical, volcano, save, gsea), `tl/` (testing, colocalization, deconvolution, io, score_signature, gene_sets, gsea), `memory/` (profiling, gpu), `qc/` (metrics, plots, spatial).
- **Gene set scoring redesign (Phase 3.5b):** Full overhaul. (1) `sc_tools/data/hallmark_human.json`: 50 bundled MSigDB Hallmark sets (offline). (2) `sc_tools.tl.gene_sets`: `load_hallmark`, `load_msigdb_json`, `load_gmt`, `list_gene_sets`, `validate_gene_signatures`, `merge_gene_signatures`, `update_gene_symbols`, `save_gene_signatures`. (3) `score_signature(method=...)`: scanpy (default), ucell (pyucell), ssgsea (gseapy). (4) `sc_tools.tl.gsea`: `run_ora` (Fisher exact + BH), `run_gsea_pseudobulk` (prerank). (5) `sc_tools.pl.gsea`: `plot_gsea_dotplot`. Optional deps: `pip install sc-tools[geneset]`. 80 tests pass.
- **Generic deconvolution module:** `sc_tools.tl.deconvolution()` with backend registry (cell2location, tangram, destvi), `extract_reference_profiles()`, per-library backed loading, memory profiling. Output: `obsm['cell_type_proportions']`, `obs['{method}_argmax']`. 14 unit tests. Thin project wrappers for ggo_visium and robin.
- **Deconvolution spatial plots:** `scripts/plot_deconvolution_spatial.py` saves per-library PNGs at 300 DPI alongside PDF. Method-specific output dirs (`figures/deconvolution/{method}/`). Spot sizing `120000/n_spots`; 98th percentile vmax; white background. Both Snakefiles updated (Phase 3.5b rules).
- **projects layout:** `visium/`, `visium_hd/`, `xenium/`, `imc/`, `cosmx_1k/`, `cosmx_6k/`, `cosmx_full_library/`; each project has `data/`, `figures/`, `metadata/`, `scripts/`, `results/`, `outputs/`.
- **create_project.sh:** `./projects/create_project.sh <project_name> <data_type>`.
- **Makefile project-aware:** `PROJECT ?= projects/visium/ggo_visium`; paths use `$(PROJECT)/...`.
- **Metadata moved:** Root `metadata/` moved to `projects/visium/ggo_visium/metadata/`.
- **Legacy merged:** `scripts/old_code/` (read-only reference).
- **Spatial multipage refactor (ggo_visium):** Logic moved from repo-root `scripts/spatial_multipage.py` into `projects/visium/ggo_visium/scripts/spatial_multipage_colocalization.py` using sc_tools: scores from obsm via `get_signature_df`, `multipage_spatial_pdf` and `truncated_similarity` in sc_tools; outputs to `figures/manuscript/macrophage_localization/` and `figures/manuscript/neutrophil_cytotoxic_tcell_localization/`. Repo-root script deprecated (docstring). Makefile target: `make spatial-multipage-colocalization`.
- **Phase 3 preprocessing module (sc_tools.pp):** Modality-aware preprocessing with GPU auto-detection (rapids-singlecell/scanpy). Normalization (`normalize_total`, `log_transform`, `scale`, `arcsinh_transform`, `filter_genes_by_pattern`, `backup_raw`). Integration (`run_scvi`, `run_harmony`, `run_cytovi`; all soft deps). Reduction/clustering (`pca`, `neighbors` with auto use_rep detection, `umap`, `leiden`, `cluster`, `run_utag`). Recipes: `preprocess(modality="visium"|"visium_hd"|"xenium"|"cosmx"|"imc")`. 38 unit tests pass, lint clean. skills.md updated with modality-specific preprocessing standards.
- **Phase 0 ingestion module (sc_tools.ingest):** Batch manifest system (`config.py`), Space Ranger / Xenium / IMC command builders (`spaceranger.py`, `xenium.py`, `imc.py`), modality-specific AnnData loaders (`loaders.py`). Snakemake Phase 0 rules and validation sentinels in robin, ggo_visium, and create_project.sh template. 26 unit tests pass, lint clean.
- **Checkpoint validation (sc_tools.validate):** `validate_p1` through `validate_p4` with auto-fix support. CLI wrapper `scripts/validate_checkpoint.py` for Snakemake integration. Validation sentinels enforce Architecture.md Section 2.2 metadata contracts. 30 unit tests pass, lint clean.
- **IMC segmentation benchmark pipeline (sc_tools.bm + sc_tools.data.imc.benchmark):** Comprehensive 4-strategy benchmark framework. Data: `sc_tools/data/imc/benchmark/` (catalog for 8 internal + 5 public datasets, BenchmarkConfig with YAML I/O, ROI validation, IMC intensity normalization, DNA channel extraction, probability map generation). Strategies: (1) Ilastik prob map + DL, (2) DNA-only + prob map + DL, (3) HuggingFace pretrained (CellViT/SAM/HoVer-Net), (4) SegFormer trainable N-channel. Infrastructure: runner with resume, SLURM job generation, CLI (`python -m sc_tools.bm.cli`, 8 subcommands), cross-dataset analysis (Wilcoxon + BH). Extended metrics: panoptic quality (PQ = SQ x DQ), boundary F1, cell type preservation. 7 new plot functions, HTML report template, shared post-processing (watershed, area filter, hole fill). `benchmark-extended` dep group in pyproject.toml. 21 new files, 8 modified, 66 new tests (394 total pass), lint clean.
- **Publication figure guidelines (skills.md §12):** Journal-specific standards for Nature/Cell/Science families. Dimension/DPI/font tables, Okabe-Ito and Paul Tol color-blind safe palettes, statistical reporting requirements, line weights, scale bars, marsilea for composite figures, pre-submission checklist. `pyproject.toml [viz]` extra for marsilea.
- **skills.md reorganization:** Restructured into 4 logical parts: (I) Analysis Pipeline §1-6, (II) Spatial and Integrative Analysis §7-9, (III) Statistical Analysis and Visualization §10-12, (IV) Engineering and CI §13-17. Colors (old §11.5) merged with palettes into §11. References consolidated into Appendix.

### In Progress — Priority order

**1. Automated cell typing + post-celltyping QC report — DONE (2026-03-06)**
- [x] `sc_tools.tl.celltype` subpackage (Plan A): ScType, gating, CellTypist, ensemble; 37 tests.
- [x] `generate_post_celltyping_report()` overhaul + tabbed navigation (Plan B): dedicated template, `qc_celltype_abundance()`, tab wrapping reports 2-4; 17 tests.
- [x] `apply_celltype_map()` consolidates project-local `apply_leiden_mapping()` into sc_tools.

**2. IMC ingestion pipeline (next)**
- [ ] Wire end-to-end IMC Snakemake workflow (Phase 0a → 0b) for lymph_dlbcl and ggo-imc.
- [x] Create IMC batch manifest templates under `metadata/phase0/`. ggo_human: `batch1_samples.tsv` (8 samples). ggo-imc: `scripts/generate_manifests.py` for auto-generating panel manifests.
- [ ] Validate against existing processed data on HPC. Both IMC projects have non-standard per-ROI h5ad layouts (not standard `cells.h5ad`).

**2. PyPI deployment (CI/CD step 4) — DONE**
- [x] Ensure `pyproject.toml` builds clean wheel/sdist.
- [x] GitHub Action for publishing on release (trusted publishing).
- [x] LICENSE file (MIT, Yoffe Lab). `requires-python >= 3.11`. Package-data for bundled JSON/HTML.
- [x] One-time PyPI setup: trusted publisher configured on pypi.org; v0.1.0 released.

**3. Registry schema extension + Phase DAG — DONE**
- [x] `sc_tools/pipeline.py`: 10-slug DAG (`PhaseSpec`, `STANDARD_PHASES`, `extend_dag`, `get_available_next`, `get_phase_checkpoint`, `validate_dag`).
- [x] Registry: `domain`/`imaging_modality` on projects; `file_role`/`validated`/`n_obs`/`n_vars` on datasets; `project_phases` table with composite PK.
- [x] Alembic migrations: `0001_initial_schema.py` (baseline), `0002_tech_taxonomy_and_phase_tracking.py`.
- [x] MCP server: `get_phase_status`, `set_phase_status` tools added.
- [x] README.md Mermaid diagram fixed (TD layout, no nested subgraphs, new phase slug names).
- [x] Architecture.md and CLAUDE.md updated with slug-to-old-code mapping.
- [x] 59 new tests (500 total pass, 12 skipped), lint clean.

**4. Phase 0b checkpoint script — PARTIALLY DONE**
- [x] Generic `scripts/ingest.py` that reads `all_samples.tsv`, calls modality loader, saves per-sample `adata.ingested.h5ad`. 5 tests.
- [ ] Snakemake rule `adata_ingested` with sentinel `data/{sample_id}/.adata.ingested.done`.

**4. Testing**
- [ ] ggo_visium project tests: `projects/visium/ggo_visium/tests/`. Integration/smoke tests.
- [ ] sc_tools package tests: `sc_tools/tests/`. Expand unit tests for pl, tl, qc.

**5. Scripts cleanup**
- [x] Project Makefiles deleted (ggo_visium, robin) — fully superseded by per-project Snakefiles.
- [x] Root Makefile retained for developer tooling only (lint, format, test, docs).
- [ ] Update scripts to use `$(PROJECT)/metadata/` paths.
- [ ] Modular scripts: config-driven; import from `sc_tools`; thin orchestration only.

### CI/CD Roadmap

Sequential; each step builds on the previous. See skills.md Sections 13 and 16.

| Step | Task | Status | Notes |
|------|------|--------|-------|
| 1 | **Linting** | Done | Ruff in `pyproject.toml`; `make lint` runs `ruff check` + `ruff format --check`. |
| 2 | **Snakemake** | Done | Per-project Snakefiles. Apptainer (Linux/HPC) / Docker (macOS) via `run_container.sh`. |
| 3 | **Sphinx docs** | Done | `docs/` with pydata-sphinx-theme + myst-nb; `make docs`; `.readthedocs.yaml`; zero warnings. |
| 4 | **PyPI deployment** | Done | `pyproject.toml` wheel/sdist; `.github/workflows/publish.yml` (trusted publishing). LICENSE added. |
| 5 | **GitHub Actions** | Done | `.github/workflows/ci.yml`: lint, test matrix (py3.10/3.11 x ubuntu/macos), docs, Snakemake dry-run, Docker build. `docs.yml` absorbed. |
| 6 | **Storage + Registry + MCP** | Done | `sc_tools/storage.py` (fsspec URI abstraction), `sc_tools/registry.py` (SQLAlchemy SQLite/PG), `sc_tools/mcp/` (FastMCP servers), `.mcp.json`. 441 tests pass. |
| 7 | **Registry schema extension + Phase DAG** | Done | `sc_tools/pipeline.py` (10-slug DAG, `extend_dag`, `get_available_next`, `get_phase_checkpoint`). Registry: `domain`/`imaging_modality` on projects; `file_role`/`validated`/`n_obs`/`n_vars` on datasets; `project_phases` table. Alembic migrations `0001`+`0002`. MCP: `get_phase_status`, `set_phase_status`. README Mermaid fixed. 500 tests pass. |

### To Do (later)

- [ ] Organize production scripts by phase within project.
- [ ] Update imports to use `sc_tools.pl.*`, `sc_tools.tl.*` everywhere.
- [ ] CosMx loader: `load_cosmx_sample()` for 1k/6k/full_library panels (deprioritized; no active CosMx project).
- [ ] Documentation: migration guide, docstrings, API docs.

### Operational notes

- **API:** `st.pl.*`, `st.tl.*`, `st.qc.*`, `st.memory.*` (scanpy-style).
- **New project:** `./projects/create_project.sh <project_name> visium|visium_hd|visium_hd_cell|xenium|imc|cosmx`. Run make with `PROJECT=projects/<type>/<name>`.
- **Legacy:** `scripts/old_code/` is read-only. Do not modify; refactor into new scripts or `sc_tools` if needed.

---

## 7. Reference

- **README.md** — Pipeline workflow diagram, phase summary.
- **Architecture.md** — Directory layout, phase details, project-specific paths, testing structure.
- **skills.md** — Coding and statistical standards.
- **Per-project:** `projects/<data_type>/<project_name>/Mission.md` and `Journal.md`.
