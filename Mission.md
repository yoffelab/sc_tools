# Mission: sc_tools Toolkit and Pipeline

**Scope:** Repository-level and generalizable goals. Project-specific objectives (e.g. TLS, macrophage analysis) live in each project's `Mission.md` under `projects/<data_type>/<project_name>/Mission.md`.

**Last Updated:** 2026-03-04 (IMC loader API fix: processed_dir; CosMx required columns)

---

## 1. Toolkit and Pipeline (General)

- **Pipeline phases (0–7):** Non-linear workflow with human-in-loop steps. See `README.md` for the diagram and phase summary.
- **sc_tools:** Reusable package (`sc_tools.pl`, `sc_tools.tl`, `sc_tools.qc`, etc.) for QC, plotting, testing, colocalization, I/O. Keep generic; project-specific logic stays in project scripts.
- **Reproducibility:** Makefile is project-aware (`PROJECT ?= projects/visium/ggo_visium`). Each project has `data/`, `figures/`, `metadata/`, `scripts/`, `results/`, `outputs/` under `projects/<platform>/<project_name>/`.
- **Standards:** All projects follow `skills.md` (statistics, significance bars, FDR). Documentation avoids apostrophes.

### Phasing scheme (sc_tools)

| Phase | Name | Notes |
|-------|------|--------|
| **0** | Upstream Raw Data Processing | (0a) Space Ranger / Xenium Ranger / IMC pipeline on HPC → `data/{sample_id}/outs/`. (0b) Load per-sample into AnnData/SpatialData → `data/{sample_id}/adata.p0.h5ad`. |
| **1** | QC and Concatenation | Load Phase 0 per-sample AnnData; per-sample QC; filter; concat via `concat_samples()`; QC reports. |
| **2** | Metadata Attachment | Human-in-loop unless `sample_metadata.csv`/`.xlsx` provided. |
| **3** | Preprocessing | Backup `adata.raw`; filter; normalize; batch correct (e.g. scVI); cluster. No automated cell typing. |
| **3.5** | Demographics | Branch from Phase 3 (parallel to 3.5b). Cohort stats, Figure 1. |
| **3.5b** | Gene scoring, automated cell typing, deconvolution | Branch from Phase 3 (parallel to 3.5). Hallmark + project signatures (target: `adata.obsm['sig:...']`); automated cell typing; optional deconvolution. Connects to Phase 4. |
| **4** | Manual Cell Typing | Human-in-loop, iterative. Skippable if automated typing (3.5b) is adequate. |
| **5** | Downstream Biology | Uses 3.5b outputs. Spatial/process analysis, figures. |
| **6–7** | Meta Analysis | Optional. Aggregate ROI/patient; downstream on aggregated. |

**Entry points:** Start at Phase 3 (preprocessed AnnData), Phase 3.5b (clustered AnnData), or Phase 4 (phenotyped AnnData). See README diagram for START HERE conditions.

**Checkpoint nomenclature:** Standard result filenames and required metadata are defined and checked in **Architecture.md** (Section 2). Use `results/adata.raw.p1.h5ad`, `adata.annotated.p2.h5ad`, `adata.normalized.p3.h5ad`, `adata.normalized.scored.p35.h5ad`, `adata.celltyped.p4.h5ad`, and `adata.{level}.{feature}.h5ad` for aggregated (level ∈ roi|patient, feature ∈ mean_expression|celltype_frequency).

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

- [x] **Modality loaders (implemented):** `sc_tools.ingest.loaders` — `load_visium_sample()`, `load_visium_hd_sample()`, `load_visium_hd_cell_sample()`, `load_xenium_sample()`, `load_imc_sample()`. Each sets `obs['sample']`, `obs['library_id']`, `obs['raw_data_dir']`, `obsm['spatial']`.
- [ ] **CosMx loader:** `load_cosmx_sample()` — reads flat CSV/Parquet files or RDS (via rpy2+anndata2ri) from NanoString/AtoMx output. Sets spatial coordinates from cell centroid (x, y in microns). Batch TSV schema: `sample_id`, `cosmx_dir`, `batch`.
- [x] **visium_hd_cell modality:** SpaceRanger 4 cell segmentation support. `load_visium_hd_cell_sample()` reads from `outs/cell_segmentation/`. Cell-level QC thresholds (like xenium). Routes to xenium recipe in `preprocess()`. `create_project.sh` accepts `visium_hd_cell` as valid data type.
- [x] **Concatenation:** `concat_samples()` merges per-sample AnnData with `calculate_qc_metrics()` applied.
- [ ] **Phase 0b checkpoint:** Per-sample `data/{sample_id}/adata.p0.h5ad` saved by `scripts/ingest.py` before concatenation. Snakemake rule `adata_p0` produces sentinel `data/{sample_id}/.adata.p0.done`.
- [ ] **SpatialData (optional):** `data/{sample_id}/spatialdata.zarr` for Visium HD and Xenium when full image pyramids / subcellular coords needed. Loader via `spatialdata-io`.

### Checkpoint Validation

- [x] **`sc_tools.validate`:** `validate_p1` through `validate_p4`, `validate_checkpoint()` dispatcher, `validate_file()` convenience. Auto-fix: `obs['batch']` to `obs['raw_data_dir']` (p1). 30 unit tests.
- [x] **CLI:** `scripts/validate_checkpoint.py` — `--phase`, `--fix`, `--warn-only`. Exit code 1 on failure.
- [x] **Snakemake validation sentinels:** `results/.adata.{name}.validated` sentinel files in robin, ggo_visium, and create_project.sh template. Downstream rules depend on sentinels.

### Phase 1: QC and Concatenation

**Input:** Per-sample `data/{sample_id}/adata.p0.h5ad` produced in Phase 0b.

- [ ] **Load Phase 0 checkpoints:** Read all per-sample `adata.p0.h5ad` files listed in `metadata/phase0/all_samples.tsv`.
- [ ] **Per-sample QC:** `sc_tools.qc.filter_spots()` (modality-aware thresholds); `compute_sample_metrics()`; `classify_samples()` (absolute + MAD outlier thresholds).
- [ ] **Concatenation:** `concat_samples()` across all passing samples → `results/adata.raw.p1.h5ad`.
- [ ] **Required annotations:** `obs['sample']`, `obs['raw_data_dir']`, `obsm['spatial']`; `X` raw counts; no normalization.

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
- [x] `generate_qc_report()`: self-contained HTML with summary cards, per-sample metrics table, embedded plots. Jinja2 template at `sc_tools/data/qc_report_template.html`.
- [x] Integrated into `scripts/run_qc_report.py` with `--modality`, `--mad-multiplier`, `--thresholds`, `--apply-filter` CLI args.
- [x] 15 new unit tests (28 total in test_qc.py); all pass, lint clean.

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
- [ ] Automated cell typing (cluster → celltype); non-transformer and transformer models.
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
- **Journal summary and Mission-as-todo workflow:** journal_summary.md at root and per project (lymph_dlbcl, ggo_visium); Mission.md is the todo list; in work mode the agent updates Mission after each prompt. Skill (`.cursor/skills/journal-and-mission-workflow/`), rule (`.cursor/rules/journal-and-mission.mdc`), Cursor settings reminder; create_project.sh creates journal_summary.md for new projects.
- **sc_tools skills as Cursor skill:** Repository root `skills.md` is exposed as Cursor skill `.cursor/skills/sc-tools-skills/SKILL.md`; agent follows skills.md for analysis and coding. Sandbox/local defaults: Docker + Snakemake (documented in skills.md §13 and in the skill).
- **sc_tools package:** `pl/` (spatial, heatmaps, statistical, volcano, save, gsea), `tl/` (testing, colocalization, deconvolution, io, score_signature, gene_sets, gsea), `memory/` (profiling, gpu), `qc/` (metrics, plots, spatial).
- **Gene set scoring redesign (Phase 3.5b):** Full overhaul. (1) `sc_tools/data/hallmark_human.json`: 50 bundled MSigDB Hallmark sets (offline). (2) `sc_tools.tl.gene_sets`: `load_hallmark`, `load_msigdb_json`, `load_gmt`, `list_gene_sets`, `validate_gene_signatures`, `merge_gene_signatures`, `update_gene_symbols`, `save_gene_signatures`. (3) `score_signature(method=...)`: scanpy (default), ucell (pyucell), ssgsea (gseapy). (4) `sc_tools.tl.gsea`: `run_ora` (Fisher exact + BH), `run_gsea_pseudobulk` (prerank). (5) `sc_tools.pl.gsea`: `plot_gsea_dotplot`. Optional deps: `pip install sc-tools[geneset]`. 80 tests pass.
- **Generic deconvolution module:** `sc_tools.tl.deconvolution()` with backend registry (cell2location, tangram, destvi), `extract_reference_profiles()`, per-library backed loading, memory profiling. Output: `obsm['cell_type_proportions']`, `obs['{method}_argmax']`. 14 unit tests. Thin project wrappers for ggo_visium and robin.
- **Deconvolution spatial plots:** `scripts/plot_deconvolution_spatial.py` saves per-library PNGs at 300 DPI alongside PDF. Method-specific output dirs (`figures/deconvolution/{method}/`). Spot sizing `120000/n_spots`; 98th percentile vmax; white background. Both Snakefiles updated (Phase 3.5b rules).
- **projects layout:** `visium/`, `visium_hd/`, `xenium/`, `imc/`, `cosmx/`; each project has `data/`, `figures/`, `metadata/`, `scripts/`, `results/`, `outputs/`.
- **create_project.sh:** `./projects/create_project.sh <project_name> <data_type>`.
- **Makefile project-aware:** `PROJECT ?= projects/visium/ggo_visium`; paths use `$(PROJECT)/...`.
- **Metadata moved:** Root `metadata/` moved to `projects/visium/ggo_visium/metadata/`.
- **Legacy merged:** `scripts/old_code/` (read-only reference).
- **Spatial multipage refactor (ggo_visium):** Logic moved from repo-root `scripts/spatial_multipage.py` into `projects/visium/ggo_visium/scripts/spatial_multipage_colocalization.py` using sc_tools: scores from obsm via `get_signature_df`, `multipage_spatial_pdf` and `truncated_similarity` in sc_tools; outputs to `figures/manuscript/macrophage_localization/` and `figures/manuscript/neutrophil_cytotoxic_tcell_localization/`. Repo-root script deprecated (docstring). Makefile target: `make spatial-multipage-colocalization`.
- **Phase 3 preprocessing module (sc_tools.pp):** Modality-aware preprocessing with GPU auto-detection (rapids-singlecell/scanpy). Normalization (`normalize_total`, `log_transform`, `scale`, `arcsinh_transform`, `filter_genes_by_pattern`, `backup_raw`). Integration (`run_scvi`, `run_harmony`, `run_cytovi`; all soft deps). Reduction/clustering (`pca`, `neighbors` with auto use_rep detection, `umap`, `leiden`, `cluster`, `run_utag`). Recipes: `preprocess(modality="visium"|"visium_hd"|"xenium"|"cosmx"|"imc")`. 38 unit tests pass, lint clean. skills.md updated with modality-specific preprocessing standards.
- **Phase 0 ingestion module (sc_tools.ingest):** Batch manifest system (`config.py`), Space Ranger / Xenium / IMC command builders (`spaceranger.py`, `xenium.py`, `imc.py`), modality-specific AnnData loaders (`loaders.py`). Snakemake Phase 0 rules and validation sentinels in robin, ggo_visium, and create_project.sh template. 26 unit tests pass, lint clean.
- **Checkpoint validation (sc_tools.validate):** `validate_p1` through `validate_p4` with auto-fix support. CLI wrapper `scripts/validate_checkpoint.py` for Snakemake integration. Validation sentinels enforce Architecture.md Section 2.2 metadata contracts. 30 unit tests pass, lint clean.
- **Publication figure guidelines (skills.md §12):** Journal-specific standards for Nature/Cell/Science families. Dimension/DPI/font tables, Okabe-Ito and Paul Tol color-blind safe palettes, statistical reporting requirements, line weights, scale bars, marsilea for composite figures, pre-submission checklist. `pyproject.toml [viz]` extra for marsilea.
- **skills.md reorganization:** Restructured into 4 logical parts: (I) Analysis Pipeline §1-6, (II) Spatial and Integrative Analysis §7-9, (III) Statistical Analysis and Visualization §10-12, (IV) Engineering and CI §13-17. Colors (old §11.5) merged with palettes into §11. References consolidated into Appendix.

### In Progress — Next immediate steps

**Testing (order: 1st ggo_visium, 2nd sc_tools, 3rd functions)**
- [ ] ggo_visium project tests: `projects/visium/ggo_visium/tests/`. Integration/smoke tests for Makefile and scripts.
- [ ] sc_tools package tests: `sc_tools/tests/`. Unit tests for pl, tl, qc with synthetic fixtures.
- [ ] Add pytest to pyproject.toml / requirements.

**Makefile & scripts**
- [ ] Align Makefile targets with Phases 0–7 (see README Pipeline Workflow).
- [ ] Update scripts to use `$(PROJECT)/metadata/` instead of `metadata/`. Affected: score_gene_signatures.py, tumor_differences.py, tls_analysis.py, manuscript_spatial_plots.py, celltyping.py, etc.
- [ ] Modular scripts: config-driven; import from `sc_tools`; thin orchestration only.

### CI/CD Roadmap

Sequential; each step builds on the previous. See skills.md Sections 13 and 16.

| Step | Task | Status | Notes |
|------|------|--------|-------|
| 1 | **Linting** | Done | Ruff in `pyproject.toml`; `make lint` runs `ruff check` + `ruff format --check`. |
| 2 | **Snakemake** | Done | Per-project Snakefiles. Apptainer (Linux/HPC) / Docker (macOS) via `run_container.sh`. |
| 3 | **Sphinx docs** | Done | `docs/` with pydata-sphinx-theme + myst-nb; `make docs`; `.readthedocs.yaml`; zero warnings. |
| 4 | **PyPI deployment** | Next | `pyproject.toml` wheel/sdist; GitHub Action on release (trusted publishing). |
| 5 | **GitHub Actions** | Pending | `.github/workflows/tests.yml`: lint + pytest + Snakemake dry-run. Last (complex setup). |

### To Do (later)

- [ ] Organize production scripts by phase within project.
- [ ] Update imports to use `sc_tools.pl.*`, `sc_tools.tl.*` everywhere.
- [ ] imc-analysis: review, identify functionalities, integrate into `sc_tools`.
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
