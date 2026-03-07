# Research Journal & Decision Log: sc_tools Repository

This journal documents **repository-level** technical and structural decisions. Project-specific analysis decisions (e.g. GGO Visium TLS, macrophage, tumor differences) are in each project's `Journal.md` under `projects/<data_type>/<project_name>/Journal.md`.

## Repo-wide notes
* **AI-assisted workflow:** Agentic scientific programming is used; project content lives under `projects/<data_type>/<project_name>/`.
* **Constraint:** Generated documentation and communications avoid apostrophes.

---

## Log Entries (toolkit / repo structure)

### [2026-03-06] - Automated cell typing (sc_tools.tl.celltype) + post-celltyping QC report overhaul

**Plan A: sc_tools.tl.celltype subpackage (automated cell typing)**

Implemented a new `sc_tools/tl/celltype/` subpackage following the same Protocol-based backend registry pattern as `sc_tools.tl.deconvolution`. Key design decisions:

- **Subpackage not flat file**: 7+ backends would exceed 700 lines. Subpackage follows `sc_tools.bm` precedent and allows deferred imports per backend.
- **Protocol + registry**: `CelltypeBackend` Protocol + `_BACKENDS` dict + `register_celltype_backend()`. Allows external code to register custom backends without modifying sc_tools.
- **Cluster-level assignment as default**: Spatial data (Visium, IMC) has noisy spots; cluster-level majority voting is the production standard. Per-cell assignment available via `assign_by="cell"` (ScType) or `majority_voting=False` (CellTypist).
- **Soft deps inside run() body**: Never at module top-level, so importing `sc_tools.tl.celltype` never fails.
- **Storage schema**: `obs[f'celltype_auto_{method}']` (Categorical), `obs[f'celltype_auto_{method}_score']` (float64), `obsm[f'celltype_proba_{method}']` (optional, DataFrame), `uns[f'celltype_auto_{method}']` (metadata dict with date, version, markers).
- **`apply_celltype_map()`**: Absorbs project-local `apply_leiden_mapping()` from robin into sc_tools. Accepts dict or JSON path, writes `obs[celltype]`, `obs[celltype_broad]`, `uns[f'{key}_colors']`.
- **Backends implemented**: ScType (pure numpy), custom hierarchical gating (DFS), CellTypist (soft dep), ensemble (majority_vote, weighted_vote). Stubs for scArches, Geneformer, scGPT, SingleR.
- **Bundled data**: `marker_db/immune_human_rna.json`, `immune_human_protein.json`, `protein_aliases.json`; `gating_templates/immune_broad.yaml`, `lymph_node.yaml`.
- **pyproject.toml**: `[celltyping]` (celltypist, scarches), `[foundation]` (transformers, datasets, accelerate) extras added.
- **Tests**: 37 new tests (TDD first); 574 total pass, 12 skipped, lint clean.

**Plan B: post-celltyping QC report overhaul + tabbed navigation**

- **New template** `post_celltyping_qc_template.html`: purple header (#6c3483) to distinguish from other reports. Sections: celltype composition (stacked bar + table), Phase 3 vs Phase 4 benchmark (bio column bold/primary), radar, UMAP grid, celltype distribution per sample, optional marker dotplot, segmentation, batch-vs-bio scatter.
- **`ranking_rows` bug fix**: `generate_post_celltyping_report()` was using `_POST_INTEGRATION_TEMPLATE` and never building `ranking_rows` — table was empty in report 4. Fixed by building ranking sorted by `bio_score` (primary for Phase 4 with validated celltypes).
- **New params**: `comparison_df_p3` (Phase 3 leiden-based metrics for side-by-side), `marker_genes` (optional dotplot dict).
- **`qc_celltype_abundance()`**: Stacked bar chart in `sc_tools/pl/qc_plots.py`; uses `adata.uns[celltype_colors]` or Okabe-Ito palette; exported from `sc_tools.pl`.
- **Tabbed navigation** (Option B — inline embed): At generation time, glob `output_dir` for prior-phase reports, extract body content + head CSS via regex, inject as tab panels. Report 2 embeds report 1; report 3 embeds 1+2; report 4 embeds 1+2+3. All CSS/JS inline — fully self-contained HTML.
- **4 new helpers** in `report_utils.py`: `_find_latest_report`, `_extract_body_content`, `_extract_head_css`, `_wrap_with_tabs`. No external deps (stdlib `re` only). Fallback: if no prior reports found, standalone HTML written as before.
- **Report 3 template improvements**: Yellow callout box for bio metrics disclaimer, `#e8f4f8` Batch column header to signal it is primary metric.
- **Snakemake**: `qc_post_celltyping` rule added to robin and ggo_visium Snakefiles; `qc_report` aggregator updated.
- **CLI**: `--marker-genes` and `--comparison-df-p3` args added to `run_qc_report.py`.
- **Tests**: 17 new tests (TDD first); 574 total pass, lint clean.

### [2026-03-06] - Option A: dual phase-name support in validate.py + CLI + Snakefile template

- **Action:** Updated three files to accept both old p-codes (deprecated, old nomenclature) and new semantic slugs simultaneously.
- **`sc_tools/validate.py`:**
  - Added `_LEGACY_PHASE_MAP = {"p1": "qc_filter", "p2": "metadata_attach", "p3": "preprocess", "p35": "scoring", "p4": "celltype_manual"}`.
  - `validate_checkpoint()` now resolves legacy codes -> slugs and emits `DeprecationWarning` (message names both old code and new slug, references Architecture.md).
  - Module docstring lists all accepted slugs and maps them to old codes explicitly.
  - All `validate_p1` through `validate_p4` docstrings marked "deprecated old nomenclature; use slug via `validate_checkpoint`".
- **`scripts/validate_checkpoint.py`:**
  - `--phase` choices expanded to all 10 values (5 slugs + 5 legacy p-codes).
  - Help text and epilog clearly show slug-to-old-code mapping and note p-codes are deprecated.
  - Explicit stderr message when legacy p-code detected in argv before delegating to validate.py.
- **`projects/create_project.sh` Snakefile template:**
  - Added phase-naming convention comment block at top of generated Snakefile.
  - All rule names updated to slug names: `adata_qc_filter`, `validate_qc_filter`, `qc_filter`, `adata_metadata_attach`, `validate_metadata_attach`, `metadata_attach`, `adata_preprocess`, `validate_preprocess`, `preprocess`, `adata_scoring`, `validate_scoring`, `scoring`, `adata_celltype_manual`, `validate_celltype_manual`, `celltype_manual`.
  - All checkpoint filenames updated to new format: `adata.raw.h5ad`, `adata.annotated.h5ad`, `adata.normalized.h5ad`, `adata.scored.h5ad`, `adata.celltyped.h5ad`.
  - All `--phase` CLI args updated to slugs.
  - Mission.md template phase table updated to slug names.
  - CLAUDE.md template Key Files table updated to new checkpoint filenames.
- **Tests:** 22 new tests in `TestValidateCheckpoint` (5 slug-acceptance tests, 2 deprecation-warning tests, 1 unknown-slug-raises test). 522 total pass.

---

### [2026-03-06] - Registry schema extension: technology taxonomy + phase DAG + per-phase tracking

- **Action:** Implemented three interlocking extensions to the registry and pipeline design: (1) `sc_tools/pipeline.py` phase DAG, (2) technology taxonomy columns on projects/datasets tables, (3) `project_phases` per-phase tracking table + Alembic migrations.
- **Phase DAG (`sc_tools/pipeline.py`):**
  - `PhaseSpec` dataclass: `label`, `depends_on`, `branch`, `checkpoint`, `optional`, `iterative`.
  - `STANDARD_PHASES` dict with 10 slugs replacing the p0a/p0b/p1/p2/p3/p35/p4 naming scheme: `ingest_raw`, `ingest_load`, `qc_filter`, `metadata_attach`, `preprocess`, `demographics`, `scoring`, `celltype_manual`, `biology`, `meta_analysis`.
  - `extend_dag(slug, spec)` for project-specific or experimental phases (no library change needed).
  - `get_available_next(completed)` respects `iterative=True` (celltype_manual can be re-entered).
  - `get_phase_checkpoint(slug, **kwargs)` returns the expected output path for a phase.
  - `validate_dag()` checks all `depends_on` references.
- **Registry schema additions (`sc_tools/registry.py`):**
  - `projects` table: added `domain` (spatial_transcriptomics, spatial_proteomics, imaging, single_cell, bulk) and `imaging_modality` (brightfield, fluorescence, multiplexed_fluorescence, probe_based, mass_spec_imaging, sequencing_based).
  - `datasets` table: added `file_role` (primary, supplementary, entry_point, spatialdata, image, metadata), `validated` (bool), `n_obs`, `n_vars`.
  - New `project_phases` table: composite PK `(project_id, phase)`, upserted on each update. Tracks status (not_started, in_progress, complete, failed, skipped), entry_phase flag, primary_dataset_id FK, n_obs/n_vars/n_samples, notes, timestamps.
  - New methods: `upsert_phase`, `get_phase`, `list_phases`, `mark_dataset_validated`. Updated `mark_phase_complete` to also upsert project_phases. Updated `status()` to return per-project phase count summary.
- **Alembic migrations:**
  - `0001_initial_schema.py`: baseline (stamp target for existing DBs; creates original 4 tables).
  - `0002_tech_taxonomy_and_phase_tracking.py`: adds new columns and `project_phases` table. Uses `batch_alter_table` for SQLite compatibility.
- **MCP server (`sc_tools/mcp/registry_server.py`):** Two new tools: `get_phase_status`, `set_phase_status`. Updated `registry_status` to include phase counts. Updated `register_dataset` with `file_role`/`validated`/`n_obs`/`n_vars` params. Updated `mark_phase_complete` to sync both legacy list and project_phases table.
- **Documentation:** Fixed README.md Mermaid diagram (TD layout, no nested subgraphs, new phase slug names). Updated phase tables in README.md, Architecture.md Section 2.1, and CLAUDE.md. Added Phase DAG reference and slug-to-old-code mapping table in Architecture.md Section 2.
- **Tests:** 59 new tests (500 total pass, 12 skipped). 33 in `test_pipeline.py`, 26 new in `test_registry.py` (3 new classes: TestProjectTaxonomy, TestDatasetFileRole, TestProjectPhases). All lint clean.
- **Key design decision:** Phase slugs stored as strings in the DB; the DAG lives in Python code. This means adding a new phase requires only adding a dict entry in `pipeline.py` or calling `extend_dag()`. No migration needed for new phases.

### [2026-03-06] - Integration benchmark (9 methods) on ggo_human IMC + 4th QC report design

- **Action:** Ran comprehensive 9-method integration benchmark on ggo_human IMC data (~248K cells, 38 markers, 7 samples, 15 celltypes) on brb HPC. Documented integration selection workflow and added 4th QC report (post-celltyping).
- **Methods benchmarked (3 normalization groups):**
  - Group A (arcsinh): Unintegrated PCA, Harmony, CytoVI (arcsinh), ComBat
  - Group B (raw counts): scVI, scANVI (semi-supervised)
  - Group C (log1p z-score): IMC Phenotyping, IMC Pheno+Harmony, Z-score+Harmony
- **Results (overall score = 0.6 * batch + 0.4 * bio):**
  1. Z-score+Harmony: 0.634 (batch 0.967, bio 0.412)
  2. IMC Pheno+Harmony: 0.633 (batch 0.968, bio 0.410)
  3. IMC Phenotyping: 0.627 (batch 0.973, bio 0.397)
  4. scANVI (raw): 0.621 (batch 0.975, bio 0.385)
  5. ComBat: 0.568 | 6. PCA: 0.564 | 7. scVI: 0.555 | 8. Harmony: 0.553 | 9. CytoVI: 0.486
- **Key findings:**
  - Z-score (log1p → per-ROI truncated z-score at ±3SD) + Harmony matches the ElementoLab phenotyping pipeline performance while being simpler (no ComBat/BBKNN).
  - scANVI with raw counts performs well (4th) despite IMC data being protein-level; its semi-supervised nature leverages celltype labels.
  - CytoVI last because scvi-tools CytoVI was not available; fell back to scVI on arcsinh-transformed data (designed for raw counts). This is a real workflow people use; the benchmark reveals it is suboptimal for this dataset.
  - Bio metrics (ARI, NMI, ASW celltype) used leiden res=1.0 labels, creating potential circularity. Batch score should be the primary selection criterion pre-celltyping.
- **Harmony PyTorch fix:** `run_harmony()` rewritten to call `harmonypy.run_harmony()` directly (bypassing scanpy wrapper). Handles Z_corr as torch tensor (`.detach().cpu().numpy()`) with shape validation for both (n_obs, n_comps) and (n_comps, n_obs) orientations. Stores `np.float32` contiguous array.
- **Report optimization:** Added `comparison_df` parameter to `generate_post_integration_report()` and `compute_integration_section()` to skip recomputation (~1.5h saved on 248K cells).
- **Design decision — 4-report QC workflow:**
  - Pre-filter (Phase 1) → Post-filter (Phase 1-2) → Post-integration (Phase 3, **batch-focused**) → Post-celltyping (Phase 4, **full bio metrics**)
  - Batch score is primary metric for integration selection. Bio metrics become meaningful only after validated celltyping.
  - Integration benchmark stores intermediates in `results/tmp/integration_test/{method}.h5ad` for re-annotation without re-running integration.
  - Selected method recorded in `results/integration_method.txt`.
- **Files changed:** `sc_tools/pp/integrate.py` (Harmony fix), `sc_tools/qc/report.py` (comparison_df), `sc_tools/qc/report_utils.py` (comparison_df, expanded _KNOWN_EMBEDDINGS), `Architecture.md` (Sections 2.3, 2.5, Phase 3, Phase 4), `Mission.md` (Phase 3 integration workflow, Phase 4 QC report).

### [2026-03-06] - IMC Cell Segmentation Benchmark Pipeline Overhaul

- **Action:** Implemented comprehensive IMC segmentation benchmark framework (`sc_tools.bm` + `sc_tools.data.imc.benchmark`). 21 new files, 8 modified files across all 7 planned phases (A through G).
- **Data infrastructure (`sc_tools/data/imc/benchmark/`):**
  - `catalog.py`: `ROIRecord`/`IMCDatasetEntry` dataclasses, `discover_rois()`, `discover_internal_datasets()`, `build_benchmark_catalog()` — scans HPC dirs for 8 internal datasets (ggo-imc, aapc, clonevo, etc.).
  - `public.py`: 5 public dataset entries (steinbock_example, jackson_2020_breast, immuncan_2024, damond_2019_pancreas, rendeiro_2021_covid) with download and standardization to steinbock convention.
  - `prepare.py`: `validate_roi_files()`, `normalize_imc_intensity()` (percentile/arcsinh/zscore/uint8), `extract_dna_channels()`, `extract_channel_by_name()`, `generate_probability_map()` (gaussian/otsu/multiscale).
  - `config.py`: `BenchmarkConfig` dataclass with YAML I/O for all benchmark parameters.
- **Four segmentation strategies:**
  - Strategy 1: Full image + Ilastik prob map + DL (extended existing `segment.py` with `run_deepcell`, Cellpose variants, `run_all_strategy1`).
  - Strategy 2: DNA-only + prob map + DL (`strategy_dna.py`).
  - Strategy 3: HuggingFace pretrained models (`strategy_hf.py` — CellViT, SAM, HoVer-Net registry with IMC-to-uint8 normalization).
  - Strategy 4: SegFormer trainable N-channel input (`strategy_vit.py` — `IMCSegmentationDataset`, `build_segformer_model`, `train_segformer`, `predict_segformer`).
- **Runner and orchestration:**
  - `runner.py`: `run_single_roi()`, `run_benchmark()` with resume support via progress.json, `aggregate_results()`.
  - `slurm.py`: sbatch script generation (CPU for strategies 1-2, GPU for 3-4), job submission, result collection.
  - `cli.py`: CLI via `python -m sc_tools.bm.cli` with 8 subcommands (catalog, download, prepare, run, run-chunk, train-vit, evaluate, slurm).
- **Extended metrics (`segmentation.py`):** Added `compute_panoptic_quality()` (PQ = SQ x DQ), `compute_boundary_metrics()` (boundary F1 at multiple tolerances), `compute_cell_type_preservation()`.
- **Cross-dataset analysis (`analysis.py`):** `compute_generalization_matrix()`, `rank_methods_per_tissue()`, `statistical_comparison()` (paired Wilcoxon + BH), `identify_failure_cases()`.
- **Extended reporting:** HTML template with 9 sections (executive summary, strategy comparison, method ranking, per-dataset heatmap, tissue analysis, generalization, statistics, runtime, failure gallery).
- **7 new plot functions (`pl/benchmarking.py`):** `plot_strategy_comparison`, `plot_dataset_heatmap`, `plot_tissue_boxplot`, `plot_generalization_radar`, `plot_runtime_comparison`, `plot_failure_gallery`, `plot_metric_correlation_scatter`.
- **Shared post-processing (`postprocess.py`):** `semantic_to_instance()` (distance transform + watershed), `filter_masks_by_area()`, `fill_holes()`.
- **Dependencies:** Added `benchmark-extended` group to `pyproject.toml` (deepcell, torch, transformers, segmentation-models-pytorch, albumentations). HPC conda env YAML.
- **Tests:** 4 new test files (catalog, strategies, runner, analysis) + extended `test_segmentation_benchmark.py`. 66 new benchmark tests, 394 total pass, 0 failures.
- **Lint:** All files pass ruff check and ruff format.

### [2026-03-04] - PyPI deployment (CI/CD step 4)

- **Action:** Implemented PyPI deployment infrastructure for sc-tools v0.1.0.
- **Changes:** (1) Created `LICENSE` (MIT, Yoffe Lab 2026). (2) Updated `pyproject.toml`: added `requires-python >= 3.10`, `[tool.setuptools.package-data]` for bundled JSON/HTML files, modernized license field to SPDX string format (setuptools deprecation fix). (3) Created `.github/workflows/publish.yml` using trusted publishing (OIDC, no API tokens) with `pypa/gh-action-pypi-publish@release/v1`.
- **Verification:** `python -m build` produces clean wheel and sdist (no warnings). Wheel contains all 7 data files (hallmark_human.json, 6 HTML templates) and LICENSE. `make lint` passes.
- **Remaining:** One-time PyPI setup needed -- add trusted publisher for `yoffelab/sc_tools` on pypi.org before first release.
- **Next CI/CD step:** GitHub Actions tests.yml (step 5).

### [2026-03-04] - IMC image loading: IMCPanelMapper, build_imc_composite, load_he_image, IMC spatial plotting

- **Action:** Implemented Phase 0b IMC image loading and H&E injection, grounded in actual ggo_human data structure exploration before writing any code.
- **Key discovery:** IMC TIFF structure is NOT individual per-channel files (as originally assumed). The ElementoLab pipeline produces a single multi-channel stack per ROI: `{roi_id}_full.tiff` (C, H, W) with companion `{roi_id}_full.csv` (index → `MarkerName(IsotopeTag)` mapping). Additional files: `_full_mask.tiff` (cell segmentation), `_full_nucmask.tiff` (nuclear mask), `_Probabilities.tiff` (ilastik). Also confirmed two CSV formats: per-ROI `*_full.csv` (two-column: index, channel name) and project-level `channel_labels.csv` (columns: channel, Target, Metal_Tag, Atom, full, ilastik).
- **Changes — `sc_tools/ingest/imc.py`:**
  - `IMCPanelMapper`: maps protein/marker names to TIFF stack channel indices. `from_full_csv()` parses `*_full.csv` (reads `MarkerName(IsotopeTag)` format, extracts protein name + isotope tag). `from_panel_csv()` parses `channel_labels.csv` (adds isotope aliases + full/ilastik flags). `set_from_var_names()` bootstraps from adata.var_names. Resolution precedence: exact protein name → full string → isotope tag → partial substring → None (warning). `build_rgb_indices()` returns TIFF stack indices for R/G/B. `get_ilastik_indices()` returns channels with ilastik=1.
  - `build_imc_composite()`: reads single `*_full.tiff` via tifffile → arcsinh(x/5) normalized (C, H, W) float32. Builds (H, W, 3) uint8 RGB composite via percentile clipping. Optional mask/probabilities loading. Returns standard squidpy spatial dict with images/hires, images/full, scalefactors (spot_diameter_fullres=1.0 for 1 px=1 um IMC), metadata/channels, rgb_channels, rgb_indices, pixel_size_um.
- **Changes — `sc_tools/ingest/loaders.py`:**
  - `load_imc_sample()`: added `load_images`, `panel_csv`, `rgb_channels`, `image_downsample` params. When `load_images=True`, looks for `tiffs/{sample_id}_full.tiff` + `{sample_id}_full.csv`; falls back to glob `*_full.tiff` with warning. Stores result in `adata.uns['spatial'][sample_id]`. All failures warn and continue (never raise).
  - `load_he_image()`: new function. Injects any TIFF (RGB/grayscale/CHW) into `adata.uns['spatial'][library_id]` for any modality. Handles grayscale-to-RGB, uint16/float-to-uint8, CHW-to-HWC, alpha-drop, downsample.
- **Changes — `sc_tools/pl/spatial.py`:** Added `plot_imc_composite()` (renders `images['hires']` RGB with channel label title; identical API to `plot_spatial_plain_he`) and `plot_imc_channel()` (renders single channel from `images['full']` by name with inferno colormap and percentile vmax).
- **Changes — `sc_tools/ingest/__init__.py`:** Exported `IMCPanelMapper`, `build_imc_composite`, `load_he_image`.
- **Tests:** 18 new tests across `TestIMCPanelMapper`, `TestBuildImcComposite`, `TestLoadImcSampleImages`, `TestLoadHeImage`. 291 total pass, 11 skipped, lint clean.
- **Architecture.md:** Added Section 2.2b with full IMC image schema (uns structure, TIFF naming conventions, CSV formats, IMCPanelMapper description). Updated Section 2.2 adata.p0.h5ad row to note optional image storage.
- **Rationale:** Loading the actual ggo_human project data first revealed the critical difference — single multi-channel TIFF vs. individual per-channel TIFFs. The `IMCPanelMapper` design handles the `MarkerName(IsotopeTag)` naming convention robustly across both `*_full.csv` and `channel_labels.csv` input formats. Storing images in the standard squidpy `uns['spatial']` format preserves `sc.pl.spatial` compatibility at no extra cost.

### [2026-03-04] - Phase 0 split into 0a/0b; per-sample AnnData/SpatialData checkpoints; CosMx loader

- **Action:** Clarified Phase 0 architecture across all documentation. Split Phase 0 into two sub-steps and added per-sample AnnData/SpatialData as the Phase 0b checkpoint. Added CosMx loader as a pending TODO.
- **Changes:**
  - Phase 0a: platform tools (Space Ranger / Xenium / IMC) → `data/{sample_id}/outs/`
  - Phase 0b: `sc_tools.ingest.loaders` → `data/{sample_id}/adata.p0.h5ad` (required) and/or `data/{sample_id}/spatialdata.zarr` (optional; Visium HD / Xenium)
  - Phase 1 redefined as QC + concatenation from Phase 0b checkpoints (was: loading + QC + concat)
  - Architecture.md: Section 2.1 new checkpoint rows for adata.p0.h5ad and spatialdata.zarr; Section 2.2 adds Phase 0 AnnData metadata contract; Phase 0 description restructured with 0a/0b subsections and modality table; Phase 1 description updated
  - CosMx: added `load_cosmx_sample()` as TODO — reads flat CSV/Parquet or RDS (rpy2+anndata2ri), centroid coords → `obsm['spatial']`
  - skills.md §1: new Phase 0 flow section with modality loader table
  - CLAUDE.md: phase table split into 0a/0b rows
  - Mission.md: Phase 0 section split into 0a/0b; Phase 1 section updated; CosMx loader TODO added
- **Rationale:** The previous docs described loading as part of Phase 1, which conflated per-sample loading with multi-sample QC and concatenation. The 0a/0b split makes responsibilities clear: 0a = HPC platform tools, 0b = portable AnnData/SpatialData per sample, 1 = QC + concat. This also enables future Snakemake rules to checkpoint per-sample AnnData before concatenation.

### [2026-03-04] - Publication figure guidelines and skills.md reorganization

- **Action:** (1) Added publication figure production guidelines to `skills.md` (now Section 12). (2) Reorganized `skills.md` into 4 logical parts with clean numbering (no more §11.5).
- **Content added (§12):** Dimension/DPI tables per journal (Nature 89/183mm, Cell 85/174mm, Science 57/121/184mm; 300+ DPI). Typography specs (Helvetica/Arial, 5-8pt, bold lowercase panels for Nature, uppercase for Cell/Science). Color-blind safe palettes (Okabe-Ito hex codes, Paul Tol, viridis/cividis; banned jet/rainbow). Statistical reporting (exact P-values, error bars, sample sizes, test names). Line weights, scale bars. Marsilea for composite figures. Supplementary standards. Pre-submission checklist.
- **Reorganization:** Part I: Analysis Pipeline (§1-6), Part II: Spatial and Integrative Analysis (§7-9, imaging modalities moved up), Part III: Statistical Analysis and Visualization (§10-12, colors/palettes consolidated into §11), Part IV: Engineering and CI (§13-17). References consolidated into Appendix.
- **Files modified:** `skills.md` (rewrite), `pyproject.toml` (added `[viz]` optional dep for marsilea), `CLAUDE.md` (Phase 0 in phase table, §12 reference), `Mission.md` (CI/CD roadmap cleaned up as table, Phase 0 in ingestion section, updated section refs).
- **Rationale:** skills.md had grown organically with awkward numbering (§11.5) and scattered related topics (colors in §11.5, palettes in §16.3, stats in §9 and §16.4). Grouping by domain makes the document navigable. CI/CD roadmap converted from numbered bold list to table for cleaner presentation.

### [2026-03-04] - Add visium_hd_cell modality and Phase 0 batch manifests for robin

- **Action:** (1) Added `visium_hd_cell` modality across sc_tools for SpaceRanger 4 cell segmentation output. New `load_visium_hd_cell_sample()` loader in `sc_tools.ingest.loaders` reads from `outs/cell_segmentation/` directory. Added to `REQUIRED_COLUMNS` (same schema as visium_hd), `_SPOT_FILTER_DEFAULTS` (cell-level thresholds matching xenium), `_SAMPLE_THRESHOLDS`, `VALID_MODALITIES` in recipes.py. Cell segmentation data routes to the xenium recipe (normalize, log1p, HVG, scale, PCA, integration) since it is single-cell resolution. (2) Updated `create_project.sh` to accept `visium_hd_cell` as valid data type. (3) Created Phase 0 batch manifests for robin: `metadata/phase0/batch1_samples.tsv` (8 samples), `metadata/phase0/batch2_samples.tsv` (8 samples), auto-generated `all_samples.tsv` (16 samples) via `collect_all_batches()`. (4) Updated robin `config.yaml` with SpaceRanger 4 paths (transcriptome, probe set).
- **Files created:** `projects/visium_hd/robin/metadata/phase0/batch1_samples.tsv`, `batch2_samples.tsv`, `all_samples.tsv`.
- **Files modified:** `sc_tools/ingest/config.py`, `sc_tools/ingest/loaders.py`, `sc_tools/ingest/__init__.py`, `sc_tools/pp/recipes.py`, `sc_tools/qc/sample_qc.py`, `projects/create_project.sh`, `projects/visium_hd/robin/config.yaml`.
- **Rationale:** SpaceRanger 4 produces cell segmentation output alongside bin-based data. Cell-level data is fundamentally different (single-cell resolution, like Xenium/CosMx) and needs separate QC thresholds and preprocessing recipe. Batch manifests enable Snakemake Phase 0 rules for robin. 248 tests pass, lint clean.

### [2026-03-04] - Phase 0 ingestion module and checkpoint validation

- **Action:** Implemented two new modules: (1) `sc_tools.validate` — checkpoint validation for phases p1 through p4, enforcing Architecture.md Section 2.2 metadata contracts. `validate_p1` checks obs[sample], obs[raw_data_dir], obsm[spatial], X >= 0. `validate_p2` adds clinical metadata check. `validate_p3` checks integration rep, cluster col, adata.raw. `validate_p35` checks signature_score/z/report in obsm/uns. `validate_p4` adds celltype/celltype_broad. Auto-fix: renames obs[batch] to obs[raw_data_dir] when fix=True. CLI wrapper: `scripts/validate_checkpoint.py`. (2) `sc_tools.ingest` — Phase 0 upstream processing: batch manifest system (per-batch TSVs under metadata/phase0/, auto-collect to all_samples.tsv), Space Ranger / Xenium / IMC command builders, modality-specific AnnData loaders (Visium, Visium HD, Xenium, IMC) with standardized obs/obsm keys, concat_samples() with QC. (3) Updated robin and ggo_visium Snakefiles with validation sentinel rules and Phase 0 spaceranger_count rules driven by manifest. (4) Updated create_project.sh Snakefile template with validation + Phase 0. 56 tests pass, lint clean.
- **Files created:** `sc_tools/validate.py`, `scripts/validate_checkpoint.py`, `sc_tools/tests/test_validate.py`, `sc_tools/ingest/__init__.py`, `sc_tools/ingest/config.py`, `sc_tools/ingest/spaceranger.py`, `sc_tools/ingest/xenium.py`, `sc_tools/ingest/imc.py`, `sc_tools/ingest/loaders.py`, `sc_tools/tests/test_ingest.py`.
- **Files modified:** `sc_tools/__init__.py`, `projects/visium_hd/robin/Snakefile`, `projects/visium/ggo_visium/Snakefile`, `projects/create_project.sh`, `Architecture.md`, `Mission.md`.
- **Rationale:** Pipeline lacked standardized upstream processing and checkpoint validation. Validation sentinels enforce metadata contracts between phases, catching issues early. Phase 0 batch manifests support HPC-based Space Ranger runs with per-batch tracking.

### [2026-03-04] - QC HTML report layout overhaul: 6-row structured layout with log10 variants

- **Action:** (1) Rewrote QC HTML report from flat auto-fit PNG grid to structured 6-row, 2-column layout with inline plot generation. Row 1: QC 2x2 pre/post histograms. Row 2: Violin pre/post. Row 3: %MT per sample. Row 4: Cross-sample comparison bar (linear + log10). Row 5: Per-sample violin grouped (linear + log10). Row 6: Scatter matrix (full width). (2) Added `log_scale` parameter to `qc_sample_comparison_bar()` and `qc_sample_violin_grouped()` with `log10(x+1)` transform and custom tick labels (0→1, 1→10, 2→100, etc.); `pct_counts_mt` stays linear in violin. (3) Added new `qc_pct_mt_per_sample()` function: per-sample %MT box plots colored by pass/fail. (4) Replaced `_collect_plot_pngs()` (file-globbing) with inline `_fig_to_base64()` generation in `report.py`. (5) Added `adata_post` parameter to `generate_qc_report()` for pre/post comparison. (6) Updated HTML template caption font from 0.8rem/#888 to 1.0rem/#555. (7) Exported `qc_pct_mt_per_sample` from `sc_tools.qc` and `sc_tools.pl`. (8) Updated `run_qc_report.py` to pass `adata_post`. All 27 tests pass, lint clean.
- **Files modified:** `sc_tools/qc/plots.py`, `sc_tools/qc/report.py`, `sc_tools/data/qc_report_template.html`, `scripts/run_qc_report.py`, `sc_tools/qc/__init__.py`, `sc_tools/pl/__init__.py`.
- **Rationale:** Previous report dumped all QC PNGs in an unstructured grid. Structured layout pairs pre/post and linear/log10 side by side for direct comparison. Log10 scale reveals sample-to-sample variation that linear scale compresses. Inline plot generation eliminates dependency on pre-saved PNG files.

### [2026-03-03] - Sphinx documentation site (CI/CD step 3)
- **Action:** Created `docs/` with pydata-sphinx-theme, myst-nb, sphinx-design, sphinx-autodoc-typehints, sphinx-copybutton. Six API pages (pp, pl, tl, qc, memory, utils) using `.. autofunction::` inline — no autosummary stub generation (avoids slow subprocess imports). Six tutorial notebooks with synthetic data, static outputs, no execution. Added `.readthedocs.yaml` for RTD hosting and `.github/workflows/docs.yml` for CI artifact. Updated `pyproject.toml` with `[docs]` extras and `[project.urls]` Documentation URL. Updated `Makefile` with `docs`, `docs-clean`, `docs-open` targets. Updated `.gitignore` with `docs/_build/` and `docs/api/generated/`.
- **Decisions:** (1) Dropped `autosummary_generate = True` — spawning subprocesses to import scanpy/squidpy/dask took 20+ min and hung. Using `.. autofunction::` directly on module pages is equivalent and fast. (2) Kept `autodoc_mock_imports = []` (empty) — mocking cupy/torch caused secondary failures in anndata 0.11.4 (`is_cupy_importable()` returns True on mock, breaking `functools.singledispatch` registration; `find_spec("torch")` raises ValueError on mock with no `__spec__`). All soft deps use try/except so empty mock list is correct. (3) Fixed `sc_tools/memory/profiling.py` top-level `import anndata as ad` → `TYPE_CHECKING` guard to prevent anndata import at autodoc time. (4) Fixed emoji in `sc_tools/memory/gpu.py` (`✅` → `[OK]`, `⚠️` → `[WARN]`) for cross-platform compatibility. Build: `build succeeded` with zero warnings.
- **Files created:** `docs/conf.py`, `docs/index.rst`, `docs/installation.rst`, `docs/api/{index,pp,pl,tl,qc,memory,utils}.rst`, `docs/tutorials/{index.rst, 01-06.ipynb}`, `docs/_static/custom.css`, `.readthedocs.yaml`, `.github/workflows/docs.yml`.
- **Files modified:** `pyproject.toml` (+[docs] extras, Documentation URL), `Makefile` (+docs targets), `.gitignore` (+docs/_build/, docs/api/generated/), `sc_tools/memory/gpu.py` (emoji fix), `sc_tools/memory/profiling.py` (TYPE_CHECKING guard), `sc_tools/pl/spatial.py` (docstring indentation fix).

### [2026-02-28] - QC report module overhaul: sample-level metrics, pass/fail classification, HTML report
- **Action:** (1) Created `sc_tools/qc/sample_qc.py`: `filter_spots()` with modality-aware defaults (visium/visium_hd/xenium/cosmx/imc), `compute_sample_metrics()` for per-sample aggregate QC (n_spots, median counts/genes, %MT stats, n_genes_detected), `classify_samples()` with absolute thresholds + adaptive MAD-based outlier detection, `save_pass_fail_lists()`, `apply_qc_filter()` (backup + spot filter + sample removal). (2) Added 3 cross-sample comparison plots to `sc_tools/qc/plots.py`: `qc_sample_comparison_bar()`, `qc_sample_violin_grouped()`, `qc_sample_scatter_matrix()` with pass/fail color coding. (3) Created `sc_tools/qc/report.py` + `sc_tools/data/qc_report_template.html` (Jinja2): self-contained HTML report with summary cards, per-sample metrics table, embedded base64 PNG plots. (4) Expanded `scripts/run_qc_report.py` with `--modality`, `--mad-multiplier`, `--thresholds`, `--spaceranger-dirs`, `--apply-filter`/`--no-filter` CLI args; post adata now optional. (5) Updated `sc_tools/qc/__init__.py` and `sc_tools/pl/__init__.py` exports. (6) Added 15 new tests (28 total in test_qc.py): filter_spots, compute_sample_metrics, classify_samples (absolute/outlier/small cohort/all pass), save_pass_fail_lists, apply_qc_filter_backup, cross-sample plots, generate_qc_report. All 26 pass (2 skipped: optional deps), lint clean.
- **Rationale:** Existing QC produced scattered individual plots with no cross-sample comparison or automated pass/fail. Inspired by 10x Space Ranger web_summary.html. Conservative thresholds minimize false drops (spatial experiments are expensive); adaptive MAD multiplier scales by cohort size (n<10: mult=5.0, n>=40: mult=2.5). Two-level filtering (spot then sample) with backup preserves original data.

### [2026-02-27] - Deconvolution spatial plots: PNG output, Snakemake integration, and figure quality improvements
- **Action:** (1) Updated `scripts/plot_deconvolution_spatial.py` to save per-library PNGs at 300 DPI alongside the PDF, with filenames `{library_id}.png` in the same folder. (2) Added deconvolution rules to both ggo_visium and robin Snakefiles under Phase 3.5b, with method-specific output directories (`figures/deconvolution/{method}/`), configurable via `config.get("deconv_method", "cell2location")`. PNG outputs tracked via `touch()` sentinel files. (3) Fixed spot sizing formula from `5000/n_spots` to `120000/n_spots` (Visium ~34, Visium HD ~2.4). (4) Improved contrast using 98th percentile for vmax instead of max. (5) Added white background to all panels.
- **Rationale:** PDFs were too large for review; per-library PNGs at 300 DPI provide quick access to individual libraries. Snakemake integration ensures deconvolution and figures are part of the reproducible pipeline. Spot sizes and contrast needed tuning for both Visium and Visium HD datasets.

### [2026-02-26] - Generic cell-type deconvolution module in sc_tools
- **Action:** Implemented `sc_tools.tl.deconvolution()` as a generic entry point with pluggable backend registry pattern. Three backends: Cell2location (default, GPU-accelerated), Tangram (OT-based), DestVI (scvi-tools). Added `extract_reference_profiles()` for memory optimization (computes mean expression per cell type, ~100x smaller than full reference; Cell2location can use this via `cell_state_df` to skip regression training). Per-library backed loading with `sc_tools.memory.profiling` integration. Output stored in `obsm['cell_type_proportions']` (DataFrame, n_spots x n_celltypes) and `obs['{method}_argmax']`. Created 14 unit tests with synthetic data and mock backend. Thin project wrappers for ggo_visium (~30 lines) and robin (~25 lines) replace the 1066-line and 282-line scripts respectively. All tests pass (31/31), lint clean.
- **Rationale:** Eliminated duplicated memory management, reference loading, and proportion extraction logic across projects. Backend registry enables adding new methods without modifying orchestration. `extract_reference_profiles` addresses ggo_visium memory issues with its large scRNA-seq reference.

### [2026-02-26] - sc_tools skills as Cursor skill; Docker + Snakemake defaults for sandbox
- **Action:** (1) Added Cursor skill `.cursor/skills/sc-tools-skills/SKILL.md` that instructs the agent to read and follow repository root `skills.md` when doing single-cell/spatial analysis, pipelines, or sandbox runs. (2) Updated `skills.md` Section 11: new subsection "Sandbox and local runs (defaults)" requiring **Docker setup and Snakemake as defaults for all sandbox runs**; clarified that Snakemake is default for sandbox/local development and Nextflow for production/CI when adopted. (3) Appendix Workflow & Infrastructure now states defaults for sandbox/local: snakemake, docker. (4) Mission.md Completed section updated to record the new skill and skills.md defaults.
- **Rationale:** Make skills.md discoverable as a Cursor skill; align agent behavior with repo standards; ensure sandbox and local pipeline runs consistently use Docker + Snakemake unless overridden.

### [2026-02-17] - Docker + conda + UV integration; project_setup.md
- **Action:** (1) Switched Dockerfile base from `python:3.10-slim` to `continuumio/miniconda3`; added conda env `sc_tools` with Python 3.10. (2) Install sc_tools via `uv pip install --python /opt/conda/envs/sc_tools/bin/python -e ".[deconvolution]"` (uv requires explicit --python for conda envs). (3) ENV PATH set to conda env bin so CMD/bash use sc_tools env. (4) Created `project_setup.md` at repo root: build, run, per-project table (ggo_visium, robin, lymph_dlbcl), local dev, verification. (5) Renamed `environment.yml` env from ggo_visium to sc_tools. (6) Added `run_docker.sh` for Robin. (7) run_docker.sh header note; README links to project_setup.md.
- **Rationale:** Conda env sc_tools inside container for HPC/consistency; single image + runtime project selection; UV for fast install; project_setup.md centralizes setup docs.

### [2026-02-17] - Signature scoring: obsm-based API (Option B single flat key)
- **Action:** (1) Added `sc_tools.tl.score_signature` to score two-level JSON signatures with scanpy `score_genes`, storing raw and z-scored matrices in `obsm['signature_score']` and `obsm['signature_score_z']` with full-path column names (e.g. Myeloid/Macrophage_Core) for traceability. Report in `uns['signature_score_report']`. (2) Extended `sc_tools.utils.signatures`: `get_signature_columns_from_obsm`, `get_signature_df`, and `get_signature_columns` now prefer obsm when present. (3) Updated `sc_tools.tl.colocalization` and `sc_tools.pl.heatmaps` to use `get_signature_df(adata)` so downstream works with obsm or obs. (4) Robin: `run_signature_scoring.py` (p3 -> p35); ggo_visium: project `score_gene_signatures.py` uses sc_tools; Makefiles point p35 to these scripts. (5) Architecture.md Section 2.2 updated for p35 obsm key names and report.
- **Rationale:** Single flat obsm key gives one place to look and full-path column names for traceability; backward compatibility via get_signature_df fallback to obs.

### [2025-02-12] - Docker/Singularity for Nextflow pipeline; auto-config local vs HPC
- **Action:** Updated skills.md and Mission.md to require containerization for the Nextflow pipeline. Use **Docker** for local runs and **Singularity/Apptainer** for HPC. Pipeline must **auto-configure** based on environment: detect local (Docker available) vs HPC (Singularity available) and select the appropriate executor. Define container image(s); publish to registry for HPC pull.
- **Rationale:** Environment parity across local and HPC; reproducibility; HPC typically provides Singularity, not Docker; local dev uses Docker.

### [2025-02-12] - Switch from Snakemake to Nextflow
- **Action:** Reverted Snakemake adoption; adopted **Nextflow** instead for the phase-dependent pipeline (Phases 1–7). Updated Mission.md, skills.md, journal_summary.md. Nextflow chosen for scalability and reproducibility.
- **Rationale:** Nextflow offers better scalability (cloud, HPC), built-in Docker/Singularity support, and a strong reproducibility story; aligns with skills.md Workflow Managers.

### [2025-02-12] - Ruff linting workflow
- **Action:** (1) Added Ruff to `pyproject.toml` dev deps and `[tool.ruff]` config (line-length 100, select E/W/F/I/B/C4/UP, ignore E501/B008/E741/E402). (2) Ran `ruff format` and `ruff check --fix` on sc_tools; fixed B007, B905, F841, C408 violations; E741 (Moran's I) and E402 (imports) ignored for scientific scripts. (3) Root `Makefile` added with `make lint` (ruff check + format --check) and `make format`; `make test` placeholder. (4) sc_tools passes lint; scripts/ has many violations, to be added incrementally.
- **Rationale:** CI/CD Roadmap item 1; ensure code compiles, passes lint, and runs before commit (skills.md).

### [2025-02-12] - Snakemake adoption and CI/CD roadmap
- **Action:** (1) Updated skills.md Section 11 to **adopt Snakemake** for the phase-dependent pipeline (Phases 1–7); existing phase workflow maps to Snakemake rules; Makefile may coexist during migration; Snakemake runs in CI (dry-run or fixtures). (2) Added Snakemake to Section 15 CI (optional step in GitHub Actions). (3) Added CI/CD Roadmap to Mission.md with ordered checklist: 1 Linting, 2 GitHub Actions, 3 Snakemake adoption, 4 Sphinx, 5 PyPI. (4) Added snakemake to Appendix (CI/CD, Linting & Docs).
- **Rationale:** Align with skills.md Workflow Managers; adopt the documented phase workflow as Snakemake for reproducibility and CI validation; provide clear implementation order.

### [2025-02-12] - journal_summary.md and Mission-as-todo workflow
- **Action:** (1) Introduced **journal_summary.md** at repo root and per project (lymph_dlbcl, ggo_visium). It holds a short summary of Journal.md so both serve as long-term memory while the summary shortens context. (2) Mission.md is the single **todo list** for progress; keep it updated with checkboxes and current status. (3) In **work mode** (executing tasks), the agent updates Mission.md after each prompt; in plan mode, optional. (4) Created skill `.cursor/skills/journal-and-mission-workflow/SKILL.md`, rule `.cursor/rules/journal-and-mission.mdc`, and added a reminder in Cursor User settings. (5) Updated Architecture.md and projects/README.md; create_project.sh now creates journal_summary.md for new projects.
- **Rationale:** Shorter context via journal summary; consistent progress tracking via Mission; explicit workflow so the agent keeps Mission current when doing concrete work.

### [2025-02-11] - ggo_visium checkpoint migration to standard names
- **Action:** Renamed checkpoint files in ggo_visium: adata.annotation.masked.h5ad → adata.annotated.p2.h5ad; adata_vg_filtered.h5ad → adata.normalized.p3.h5ad; adata.img.genescores.h5ad → adata.normalized.scored.p35.h5ad. Backed up originals to results/backup260211/. Updated all scripts (annotation2mask_img, preprocessing, score_gene_signatures, signature_heatmap*, manuscript_spatial_plots, tumor_differences, macrophage_localization, process_colocalization, tls_analysis, spatial_analysis, create_tls_anndata, spatial_multipage, celltype_deconvolution_phase3, test_cell2location_spots, basic_analysis, celltype_deconvolution_tangram) and Makefile to use new paths.
- **Rationale:** Align ggo_visium with standard checkpoint nomenclature (Architecture.md Section 2).

### [2025-02-09] - Standard checkpoint nomenclature and required metadata (Architecture.md)
- **Action:** Defined mandatory checkpoint filenames and required metadata in Architecture.md Section 2. Phase outputs: `adata.raw.p1.h5ad`, `adata.annotated.p2.h5ad`, `adata.normalized.p3.h5ad`, `adata.normalized.scored.p35.h5ad`, `adata.celltyped.p4.h5ad`; aggregated (Phases 5–7): `adata.{level}.{feature}.h5ad` with level ∈ {roi, patient}, feature ∈ {mean_expression, celltype_frequency}. Section 2.2 lists required adata contents per checkpoint for validation. Operational rule added: new pipelines must use these names; validators should check metadata.
- **Rationale:** Consistent naming and metadata expectations allow tooling to validate checkpoints and scripts to assume standard paths. Legacy names allowed during migration.

### [2025-02-09] - Phase 3.5b parallel to 3.5; obsm storage TODO; diagram and Makefile
- **Action:** (1) Phase 3.5b is a separate branch from Phase 3 (parallel to 3.5), not after 3.5; it connects to Phase 4; Phase 4 is skippable if automated cell typing adequate. (2) Phase 3 no longer includes automated cell typing; that moved to Phase 3.5b. (3) Gene signature storage TODO: use `adata.obsm['sig:hallmark']` and `adata.obsm['sig:{signature_name}']` for `metadata/{signature_name}.json`; always apply Hallmark plus project-provided signatures (Mission.md, Architecture.md). (4) README diagram: 3.5 and 3.5b both branch from P3; START HERE annotated with start conditions (preprocessed / clustered / phenotyped AnnData); Phase 5 and 6–7 stacked vertically in one subgraph to save horizontal space. (5) Makefile: celltyping moved from Phase II to Phase 3.5b; phase2 now stops at scvi.leiden.h5ad; phase3.5b builds phenotyped + genescores (+ optional deconvolution).

### [2025-02-09] - Gene scoring and deconvolution moved to Phase 3.5b
- **Action:** Renamed Phase III (gene signature scoring + cell type deconvolution) to **Phase 3.5b** in phasing and phase diagram. Flow is now: Phase 3 (Preprocessing) → Phase 3.5 (Demographics) → Phase 3.5b (Gene Scoring & Deconvolution) → Phase 4 → Phase 5.
- **Files:** Architecture.md (new Phase 3.5b section; Phase 5 updated to reference 3.5b outputs), README.md (Mermaid diagram + phase table), projects/visium/ggo_visium/Makefile (phase3.5b target; phase3 kept as alias for backward compatibility).
- **Rationale:** Aligns pipeline documentation with phase numbering (3.5 = Demographics, 3.5b = Gene scoring/deconvolution before Manual Cell Typing and Downstream Biology).

### [2025-02-09] - Signature heatmaps: generic code in sc_tools, versioned figures
- **Action:** Added `scripts/signature_heatmap_versioned.py` that uses `sc_tools.pl.heatmaps.signature_score_heatmap` and `sc_tools.pl.save.save_figure` for heatmap/clustermap generation with datetime-stamped output. Fixed `sc_tools.pl.heatmaps` DataFrame construction (single dict for annotations, not generator of dicts).
- **Output:** Figures go to `figures/manuscript/signature_heatmaps/pdf/` and `.../png/` with filenames `YYDDMM.hh.mm.<basename>.pdf` (and .png) for versioning.
- **Makefile:** `phase5-direct` now runs `signature_heatmap_versioned.py`; added target `signature_heatmaps_versioned` for heatmaps only.
- **Rationale:** Minimal project-specific code in ggo_visium; generic heatmap/clustermap logic lives in sc_tools.pl for maintainability; versioned filenames improve figure tracking.

### [2025-02-09] - Pipeline Phasing Overhaul (Phases 1–7)
- **Action:** Rewrote pipeline into 7 phases with non-linear workflow and human-in-loop steps.
- **Phases:** (1) Data Ingestion & QC; (2) Metadata Attachment (HIL unless sample_metadata.csv); (3) Preprocessing; (3.5) Demographics; (4) Manual Cell Typing (HIL, iterative); (5) Downstream Biology; (6–7) Meta Analysis (optional).
- **Additions:** Created `WORKFLOW.md` with Mermaid diagram; `sc_tools.qc` placeholder (metrics, spatial, plots); project `metadata/sample_metadata.csv`, `metadata/celltype_map.json` as bypass files.

### [2025-02-09] - Merged REORGANIZATION_STATUS and REORGANIZATION_PLAN into Mission.md
- **Action:** Consolidated REORGANIZATION_STATUS.md and REORGANIZATION_PLAN.md into Mission.md Section 6 (Reorganization & Implementation Status).
- **Result:** Single source of truth for completed, in-progress, and to-do tasks. Deleted REORGANIZATION_STATUS.md and REORGANIZATION_PLAN.md.

### [2025-02-09] - Testing Strategy and Implementation Order
- **Action:** Documented two-layer testing (package + project) in Mission.md and Architecture.md.
- **Layers:** `sc_tools/tests/` (unit tests for package); `projects/<platform>/<project>/tests/` (integration tests for pipeline).
- **Implementation order:** (1) ggo_visium unit tests first, (2) sc_tools unit tests, (3) implement functions. All new code must compile and pass tests.
- **Rationale:** Establish baseline with project tests before refactoring to sc_tools; ensure guaranteed behavior of sc_tools before use in scripts.
- **GGO Visium Mission:** Unit tests added as next priority; roadmap updated.

### [2025-02-09] - Project-Specific Paths and Metadata Move
- **Action:** Moved root `metadata/` (gene_signatures.json, obs.csv) to `projects/visium/ggo_visium/metadata/`. Removed root `./metadata` folder.
- **Rationale:** All metadata, results, and figures are project-specific. No root-level metadata, results, or figures to avoid confusion with sc_tools (generic package).
- **Mission.md:** Expanded Section 2 as comprehensive to-do list with implementation status per phase.
- **Output paths:** `figures/QC/raw/`, `figures/QC/post/`; required `adata.obs['sample']`, `adata.obs['raw_data_dir']`.
- **Files updated:** Mission.md, Architecture.md, projects/visium/ggo_visium/Mission.md, projects/README.md, README.md.

### [2025-12-30] - Codebase Distillation & Archive
- **Action:** Audited scripts and moved a large portion of inherited files to a legacy folder (later merged into `scripts/old_code/` per project).
- **Rationale:** Legacy scripts were creating namespace noise; isolating them let the pipeline focus on the scVI-based workflow.
- **Decision:** Plan-Act-Reflect protocol for new analysis steps.

### [2025-12-30] - Architecture & Mission Alignment
- **Action:** Created `Architecture.md` and root `Mission.md` to guide pipeline and script roles.
- **Decision:** Legacy scripts in each project's `scripts/old_code/` are read-only references.

### [2025-02-08] - Architecture Update and Script Sanity Check
- **Action:** Updated `Architecture.md` with current layout, data flow, and script sanity check (active vs legacy vs supplementary).
- **Result:** Single place for pipeline order, script roles, and refactor rules.

### [2025-02-08] - Legacy Folder Merge: Single `scripts/old_code/`
- **Action:** Merged former `old code/` into `old_code/` per project; removed redundant folder and updated docs.
- **Result:** One read-only legacy location per project: `scripts/old_code/`.

### [2025-02-08] - Scalable Project Structure: projects/ by Data Type
- **Action:** Introduced `projects/` with data-type folders (visium, visium_hd, xenium, imc); each project has data/, figures/, metadata/, scripts/, results/, outputs/.
- **Scripts:** `create_project.sh <project_name> <data_type>`, `migrate_to_ggo_visium.sh`; Makefile project-aware (`PROJECT ?= projects/visium/ggo_visium`).
- **Result:** Repository ready for multiple projects and platforms.

### [2025-02-08] - Mission and Journal Split: Toolkit vs Project-Specific
- **Action:** Split Mission and Journal into repository-level vs per-project.
- **Rationale:** Toolkit goals (pipeline phases, sc_tools, general analysis) apply to all projects; study-specific goals (TLS, macrophage, lung tumor evolution) apply only to the relevant project.
- **Changes:**
  - **Root `Mission.md`:** Toolkit and pipeline (general); points to each project's Mission for study-specific aims.
  - **Root `Journal.md`:** Repo structure, architecture, script sanity check, legacy merge, scalable layout, and this split. No analysis-specific implementation details.
  - **Per-project:** Each project has its own `Mission.md` and `Journal.md` (e.g. `projects/visium/ggo_visium/`). Created by `create_project.sh` for new projects; GGO Visium content moved into `projects/visium/ggo_visium/Mission.md` and `Journal.md`.
- **Result:** Clear separation: root = toolkit and structure; projects = study-specific missions and analysis decision log.

### [2025-02-08] - Removed ggo_analysis and ggo_tools
- **Action:** Deleted `ggo_analysis/` and `ggo_tools/` folders.
- **Rationale:** ggo_analysis was only a placeholder (no code). ggo_tools contained only `memory/profiling.py`, which is already present in `sc_tools/memory/profiling.py` (same API: get_memory_usage, log_memory, aggressive_cleanup, estimate_adata_memory, check_memory_threshold). No merge needed; use `sc_tools.memory` for all memory utilities.
- **Result:** Single toolkit package: `sc_tools` only. Architecture.md, README.md, and .gitignore updated to drop references to ggo_tools and ggo_analysis.

### [2026-02-27] - Generic deconvolution module: NaN fix and method-specific output naming

- **Action:** Fixed `sc_tools.tl.deconvolution()` to handle already-normalised scRNA-seq references (Seurat SCTransform data with negative values). Added method-specific output file naming (`adata.deconvolution.{method}.h5ad`).
- **Root cause:** ggo_visium reference (`seurat_object.h5ad`) has negative values from SCTransform scaling (min=-11.3). The `max_val > 100` check triggered double normalisation, and `log1p` on negative values produced NaN throughout the expression matrix, causing Tangram mapping to return all-NaN.
- **Fixes in `sc_tools/tl/deconvolution.py`:**
  1. Detect negative values in reference X and skip normalisation (already scaled).
  2. Filter zero-count cells before normalisation (safety for raw count references).
  3. Replace any remaining NaN/Inf in X after normalisation (safety net).
  4. In TangramBackend: check for NaN in inputs before mapping, detect all-NaN mapping matrix and return None early.
- **Output naming:** Both wrapper scripts and Snakefiles now use `adata.deconvolution.{method}.h5ad` (e.g. `adata.deconvolution.tangram.h5ad`). Snakefiles accept `--config deconv_method=tangram` to select method.
- **Results:** ggo_visium Tangram: 29,952 spots x 31 cell types (valid, no NaN). Robin Tangram: 567,456 spots x 39 cell types (already worked). Both spatial proportion PDFs generated.
- **Note:** Cell2location is user-preferred method but requires NVIDIA GPU (not available on macOS ARM). Tangram used as CPU fallback. Wrapper defaults set to cell2location for GPU environments.

## 2026-03-06: Storage abstraction, registry DB, and MCP servers

**Goal:** Add three interlocking layers to make sc_tools work across multiple storage
platforms (HPC, S3, Box) and surface capabilities to Claude Code as callable tools.

**Storage layer (`sc_tools/storage.py`):**
- `resolve_fs(uri)`: returns `(AbstractFileSystem, path)` for any URI scheme using fsspec
- `open_file(uri, mode)`: context manager for any URI
- `with_local_copy(uri)`: downloads remote files to tmp, yields local Path (for scanpy/tifffile compat)
- `smart_read_h5ad(uri)`: reads .h5ad from local or remote URI
- `smart_write_checkpoint(adata, uri, fmt)`: writes h5ad or zarr to any URI
- `smart_read_csv(uri)`: reads CSV/TSV from any URI
- Soft deps: fsspec always available (via dask), s3fs/sshfs/gcsfs/adlfs/boxfs optional

**Registry database (`sc_tools/registry.py`):**
- SQLAlchemy ORM with 4 tables: projects, datasets, slurm_jobs, agent_tasks
- SQLite default (`~/.sc_tools/registry.db`); PostgreSQL via `SC_TOOLS_REGISTRY_URL`
- `Registry` class with full CRUD API for all tables
- CLI: `python -m sc_tools registry status`
- Migrations directory `sc_tools/migrations/` with Alembic env for future schema changes

**MCP servers (`sc_tools/mcp/`):**
- `tools_server.py`: FastMCP sc-tools server with 8 analysis tools
- `registry_server.py`: FastMCP sc-registry server with 7 bookkeeping tools
- `.mcp.json` at repo root configures both servers for Claude Code auto-discovery
- Tools in tools_server generate commands for long-running tasks (avoids blocking)

**Ingest updates:**
- `load_imc_sample()` now accepts direct `cells.h5ad` URIs (local or remote)
- `load_batch_manifest()` accepts remote URIs via `smart_read_csv()`
- `collect_all_batches()` supports remote directory globbing via fsspec

**pyproject.toml:** Added `[storage]`, `[storage-box]`, `[registry]`, `[mcp]` extras

**Results:** 441 tests pass (47 new), lint clean.
