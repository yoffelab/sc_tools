# Codebase Structure

**Analysis Date:** 2026-03-20

## Directory Layout

```
sc_tools/                          # Python package root
├── __init__.py                    # Package initialization; imports main submodules
├── __main__.py                    # CLI entry point (python -m sc_tools)
├── biodata.py                     # PlatformSpec registry for known data platforms
├── storage.py                     # fsspec URI abstraction for multi-backend storage
├── validate.py                    # Checkpoint validation (required obs/obsm per phase)
├── pipeline.py                    # Phase DAG definition and helpers
├── registry.py                    # SQLAlchemy ORM for project/phase tracking (optional)
│
├── pp/                            # Preprocessing & integration
│   ├── __init__.py               # Public API exports
│   ├── recipes.py                # Main entry point: preprocess() dispatcher
│   ├── normalize.py              # backup_raw, normalize_total, log_transform, scale
│   ├── integrate.py              # run_scvi, run_harmony, run_cytovi, run_combat, etc.
│   ├── reduce.py                 # pca, neighbors, umap, leiden, cluster, run_utag
│   ├── strategy.py               # SmallStrategy, LargeStrategy, select_strategy
│   ├── projection.py             # scVI cell projection onto reference atlas
│   ├── _gpu.py                   # GPU detection and rapids-singlecell fallback
│   └── integration_configs.py     # get_scanvi_config, get_resolvi_ss_config
│
├── qc/                            # Quality control
│   ├── __init__.py               # Public API exports
│   ├── metrics.py                # calculate_qc_metrics, filter_cells, filter_genes
│   ├── spatial.py                # spatially_variable_genes, spatially_variable_genes_per_library
│   ├── plots.py                  # qc_2x2_grid, qc_spatial_multipage, qc_violin_metrics, etc.
│   ├── sample_qc.py              # compute_sample_metrics, classify_samples, filter_by_qc
│   ├── report.py                 # generate_qc_report, generate_pre_filter_report, etc.
│   ├── report_utils.py           # Helper functions for report generation (fig_to_base64, etc.)
│   └── doublet.py                # score_doublets_solo (optional dependency)
│
├── tl/                            # Analysis tools
│   ├── __init__.py               # Public API exports
│   ├── testing.py                # Statistical testing (Mann-Whitney, FDR correction)
│   ├── colocalization.py         # Spatial colocalization (correlation, Moran's I)
│   ├── score_signature.py        # score_signature for gene set scoring
│   ├── gene_sets.py              # load_hallmark, load_msigdb_json, load_gmt, etc.
│   ├── gsea.py                   # run_ora, run_gsea_pseudobulk
│   ├── deconvolution.py          # deconvolution (cell2location, tangram, destvi)
│   ├── io.py                     # io_save (tl.io namespace)
│   └── celltype/                 # Cell type annotation
│       ├── __init__.py          # Public API exports
│       ├── annotate.py          # annotate_celltypes() dispatcher
│       ├── apply.py             # apply_celltype_map()
│       ├── _base.py             # BaseAnnotator abstract class
│       ├── _sctype.py           # ScTypeAnnotator
│       ├── _celltypist.py       # CelltypistAnnotator
│       ├── _scarches.py         # scArchesAnnotator
│       ├── _geneformer.py       # GeneformerAnnotator
│       ├── _scgpt.py            # scGPTAnnotator
│       ├── _custom_gates.py     # CustomGatesAnnotator
│       └── _ensemble.py         # EnsembleAnnotator
│
├── pl/                            # Plotting
│   ├── __init__.py               # Public API exports
│   ├── spatial.py                # Spatial plot wrappers
│   ├── heatmaps.py               # Heatmap/clustermap utilities
│   ├── statistical.py            # Statistical annotations (bars, asterisks)
│   ├── volcano.py                # Volcano plot utilities
│   ├── gsea.py                   # plot_gsea_dotplot
│   ├── qc_plots.py               # qc_celltype_abundance, qc_cluster_distribution, qc_umap_grid
│   ├── benchmarking.py           # plot_segmentation_comparison_table, plot_integration_radar, etc.
│   ├── save.py                   # save_figure
│   └── (QC plot re-exports)      # from sc_tools.qc.plots import *
│
├── ingest/                        # Phase 0: Data ingestion
│   ├── __init__.py               # Public API exports
│   ├── config.py                 # load_batch_manifest, validate_manifest, collect_all_batches
│   ├── loaders.py                # load_visium_sample, load_xenium_sample, load_imc_sample, etc.
│   ├── spaceranger.py            # build_spaceranger_count_cmd, build_batch_commands
│   ├── xenium.py                 # Xenium Ranger commands and loaders
│   ├── imc.py                    # IMCPanelMapper, build_imc_composite
│   └── slurm.py                  # build_sbatch_header, write_sbatch_script, generate_phase0_inventory
│
├── bm/                            # Benchmarking
│   ├── __init__.py               # Public API exports
│   ├── segmentation.py           # compute_morphology_metrics, score_segmentation, compare_segmentations
│   ├── integration.py            # compute_integration_metrics, compare_integrations, run_integration_benchmark
│   ├── mask_io.py                # load_mask, load_cellpose_mask, load_stardist_mask, etc.
│   ├── segment.py                # run_cellpose, run_stardist, run_deepcell, run_all_strategy1
│   ├── report.py                 # HTML report generation for benchmarking
│   ├── runner.py                 # Benchmark orchestration across ROIs
│   ├── analysis.py               # Cross-dataset statistics
│   ├── slurm.py                  # SLURM job generation
│   ├── cli.py                    # Command-line interface for benchmarking
│   ├── deepcell_runner.py        # DeepCell/Mesmer wrapper
│   ├── strategy_dna.py           # DNA-only segmentation strategy
│   ├── strategy_hf.py            # HuggingFace pretrained model strategy
│   ├── strategy_vit.py           # SegFormer trainable strategy
│   └── postprocess.py            # Shared post-processing utilities
│
├── memory/                        # Memory and GPU utilities
│   ├── __init__.py               # Public API exports
│   ├── env.py                    # Environment detection (GPU, SLURM, HPC)
│   ├── gpu.py                    # GPU memory management
│   ├── profiling.py              # Memory profiling utilities
│   └── (GPU-aware operations)    # Called by pp module
│
├── data/                          # Data resources
│   ├── __init__.py               # Public API exports
│   ├── io.py                     # load_marker_db, etc.
│   ├── *.json                    # Gene marker databases, gating templates
│   └── imc/benchmark/            # IMC benchmark data (configs, catalogs, public datasets)
│       ├── config.py
│       ├── catalog.py
│       ├── prepare.py
│       └── public.py
│
├── assets/                        # HTML/YAML templates
│   ├── qc_report_template.html
│   ├── pre_filter_qc_template.html
│   ├── post_filter_qc_template.html
│   ├── post_integration_qc_template.html
│   ├── post_celltyping_qc_template.html
│   ├── integration_report_template.html
│   └── *.yaml                    # YAML configuration templates
│
├── utils/                         # Utility functions
│   ├── __init__.py               # Public API exports
│   ├── signatures.py             # Gene signature curation helpers
│   ├── save.py                   # save_figure helpers
│   └── checkpoint.py             # Checkpoint save/load helpers
│
├── mcp/                           # Model Context Protocol (Claude integration)
│   ├── __init__.py               # Exports
│   ├── registry_server.py        # MCP tools for registry queries
│   └── tools_server.py           # MCP tools for execution
│
├── migrations/                    # Database migrations (optional, SQLAlchemy/Alembic)
│   ├── env.py
│   ├── alembic.ini
│   └── versions/                 # Migration scripts (0001_initial_schema.py, etc.)
│
├── tests/                         # Unit and integration tests
│   ├── __init__.py
│   ├── test_pp.py                # Preprocessing tests (synthetic data)
│   ├── test_pp_real_data.py      # Preprocessing tests (real project data)
│   ├── test_qc.py                # QC tests (synthetic data)
│   ├── test_qc_real_data.py      # QC tests (real project data)
│   ├── test_qc_doublet.py        # Doublet detection tests
│   ├── test_tl_real_data.py      # Tool tests (real data)
│   ├── test_gr_real_data.py      # Gene regulation tests (real data)
│   ├── test_ingest.py            # Ingestion tests
│   ├── test_ingest_cosmx.py      # CosMx-specific ingestion tests
│   ├── test_ingest_script.py     # Script-level ingestion tests
│   ├── test_deconvolution.py     # Deconvolution tests
│   ├── test_gene_sets.py         # Gene set loading tests
│   ├── test_gsea.py              # GSEA tests
│   ├── test_celltype_auto.py     # Automated cell typing tests
│   ├── test_benchmark_runner.py  # Benchmark orchestration tests
│   ├── test_benchmark_catalog.py # Benchmark catalog tests
│   ├── test_benchmark_analysis.py # Benchmark analysis tests
│   ├── test_benchmark_strategies.py # Segmentation strategy tests
│   ├── test_segmentation_benchmark.py # Segmentation benchmark tests
│   ├── test_bm_report.py         # Benchmark report generation tests
│   ├── test_integration_benchmark.py # Integration benchmark tests
│   ├── test_biodata.py           # Platform registry tests
│   ├── test_mcp_biodata.py       # MCP biodata server tests
│   ├── test_mcp_run_full_phase.py # Full phase execution tests
│   ├── test_migrations.py        # Database migration tests
│   ├── test_registry.py          # Registry ORM tests
│   ├── test_storage.py           # Storage backend tests
│   ├── test_validate.py          # Checkpoint validation tests
│   ├── test_score_signature.py   # Gene signature scoring tests
│   └── conftest.py               # pytest fixtures
│
└── gr/                            # Gene regulation (optional, advanced)
    ├── __init__.py               # Public API exports
    └── (functions)               # Functions for regulatory analysis
```

## Directory Purposes

**`sc_tools/pp/`:**
- Purpose: Preprocessing and integration orchestration
- Contains: Normalization, batch correction, dimensionality reduction, clustering, scaling strategies
- Key files: `recipes.py` (main entry), `integrate.py` (batch methods), `reduce.py` (DR/clustering), `_gpu.py` (resource detection)

**`sc_tools/qc/`:**
- Purpose: Quality control metrics, filtering, reporting
- Contains: Per-cell/per-gene metrics, spatial analysis, QC sample classification, HTML report generation
- Key files: `metrics.py` (QC calculations), `sample_qc.py` (per-sample classification), `report.py` (HTML generation)

**`sc_tools/tl/`:**
- Purpose: Analysis tools (statistical testing, deconvolution, cell typing, enrichment)
- Contains: Gene set loading, signature scoring, GSEA, colocalization, cell type annotation dispatcher
- Key files: `celltype/annotate.py` (cell typing dispatcher), `deconvolution.py` (deconvolution methods), `gene_sets.py` (gene DB loaders)

**`sc_tools/pl/`:**
- Purpose: Visualization wrappers and report figures
- Contains: Spatial plots, heatmaps, statistical annotations, benchmark plots, re-exports of QC plots
- Key files: `benchmarking.py` (comparison figures), `qc_plots.py` (post-integration plots)

**`sc_tools/ingest/`:**
- Purpose: Phase 0 ingestion and command generation
- Contains: Platform-specific loaders, batch manifest parsing, SLURM script generation
- Key files: `loaders.py` (modality-specific AnnData loaders), `config.py` (manifest parsing), `slurm.py` (HPC submission)

**`sc_tools/bm/`:**
- Purpose: Benchmarking framework for segmentation and integration comparison
- Contains: Segmentation scorers, integration metrics, model runners, HTML reporting
- Key files: `segmentation.py` (scoring logic), `integration.py` (batch correction metrics), `report.py` (figures)

**`sc_tools/data/`:**
- Purpose: Bundled biological resources (gene markers, gating templates, benchmark datasets)
- Contains: JSON gene marker databases, YAML gating rule templates, IMC benchmark catalogs
- Key files: `io.py` (loader functions), `imc/benchmark/` (IMC-specific resources)

**`sc_tools/assets/`:**
- Purpose: HTML/CSS/YAML templates for report generation
- Contains: Jinja2-rendered templates for QC reports, integration reports, post-celltyping reports
- Key files: `*_template.html` (one per report type)

**`sc_tools/utils/`:**
- Purpose: Generic helper functions
- Contains: Gene signature curation, file saving, checkpoint I/O helpers
- Key files: `signatures.py` (gene set helpers), `checkpoint.py` (save/load wrappers)

**`sc_tools/memory/`:**
- Purpose: Hardware resource detection and management
- Contains: GPU availability detection, SLURM environment parsing, memory profiling
- Key files: `env.py` (HPC environment detection), `gpu.py` (GPU memory utilities)

**`sc_tools/mcp/`:**
- Purpose: Model Context Protocol integration for Claude orchestrator
- Contains: MCP tool definitions for phase status queries, registry updates
- Key files: `registry_server.py` (registry queries), `tools_server.py` (execution tools)

**`sc_tools/migrations/`:**
- Purpose: Database schema versioning (optional, only if SQLAlchemy installed)
- Contains: Alembic migration scripts for registry database schema updates
- Key files: `versions/*.py` (one per schema change)

**`sc_tools/tests/`:**
- Purpose: Unit and integration tests
- Contains: Synthetic fixture-based tests and optional real-data tests (skipped if data unavailable)
- Key files: `conftest.py` (pytest fixtures), `test_pp.py` (preprocessing), `test_qc.py` (QC metrics)

## Key File Locations

**Entry Points:**
- `sc_tools/__init__.py`: Package initialization; imports submodules (pp, pl, tl, qc, etc.)
- `sc_tools/__main__.py`: CLI entry point for `python -m sc_tools registry status`
- `sc_tools/registry.py`: Registry ORM initialization (optional, if SQLAlchemy installed)

**Configuration:**
- `pyproject.toml` (repo root): Package metadata, dependencies, tool configuration (ruff, pytest)
- `sc_tools/pipeline.py`: Phase DAG definition (PhaseSpec instances for each phase)
- `sc_tools/biodata.py`: Platform registry (KNOWN_PLATFORMS dict)

**Core Logic:**
- `sc_tools/pp/recipes.py`: Main preprocessing dispatcher (`preprocess()` function)
- `sc_tools/qc/report.py`: QC report generation (`generate_qc_report`, `generate_pre_filter_report`, etc.)
- `sc_tools/tl/celltype/annotate.py`: Cell type annotation dispatcher (`annotate_celltypes()`)
- `sc_tools/ingest/loaders.py`: Modality-specific loaders (`load_visium_sample`, `load_xenium_sample`, etc.)

**Data Access:**
- `sc_tools/storage.py`: fsspec URI abstraction (`resolve_fs`, `open_file`, `read_h5ad`, `write_h5ad`)
- `sc_tools/validate.py`: Checkpoint validation (`validate_checkpoint`, `validate_file`)
- `sc_tools/pipeline.py`: Phase lookup (`get_phase_checkpoint`, `get_available_next`)

**Testing:**
- `sc_tools/tests/test_pp.py`: Preprocessing tests (synthetic data)
- `sc_tools/tests/test_qc.py`: QC tests (synthetic data)
- `sc_tools/tests/conftest.py`: pytest fixtures (adata_counts, adata_imc, etc.)

## Naming Conventions

**Files:**
- Implementation: `module.py` (single noun or underscore-separated compound, lowercase)
  - Examples: `recipes.py`, `normalize.py`, `report_utils.py`, `_gpu.py`
- Tests: `test_<module>.py` or `test_<component>_<scenario>.py`
  - Examples: `test_pp.py`, `test_ingest_cosmx.py`, `test_gr_real_data.py`
- Templates: `<stage>_template.html` or `<stage>_template.yaml`
  - Examples: `pre_filter_qc_template.html`, `post_integration_qc_template.html`

**Directories:**
- Module packages: lowercase singular/plural nouns
  - Examples: `pp`, `qc`, `tl`, `pl`, `ingest`, `bm`, `utils`, `memory`, `data`, `mcp`, `migrations`
- Subpackage directories: descriptive plural or functional names
  - Examples: `celltype/` (sub-package under tl), `benchmark/` (sub-package under data/imc)

**Functions:**
- Dispatch/recipe functions: verb form or action phrase, lowercase with underscores
  - Examples: `preprocess()`, `annotate_celltypes()`, `run_scvi()`, `calculate_qc_metrics()`
- Getter functions: `get_<object>` or `<object>_for_<context>`
  - Examples: `get_platform()`, `get_phase_checkpoint()`, `platform_for_project()`
- List/loader functions: `load_<resource>` or `list_<objects>`
  - Examples: `load_hallmark()`, `load_msigdb_json()`, `list_gene_sets()`, `list_platforms()`
- Boolean check functions: `is_<condition>`
  - Examples: `_is_local()` (private helper)
- Private/internal functions: prefixed with `_`
  - Examples: `_local_path()`, `_fig_to_base64()`, `_wrap_with_tabs()`

**Classes:**
- CapitalCase (PascalCase)
  - Examples: `PlatformSpec`, `PhaseSpec`, `BaseAnnotator`, `ScaleStrategy`, `SmallStrategy`
- ORM models: CapitalCase
  - Examples: `Project`, `Sample`, `Phase`, `Checkpoint` (in registry.py)

**Variables:**
- Local/module variables: lowercase with underscores
  - Examples: `adata`, `results`, `n_obs`, `n_vars`, `gene_names`
- Constants: UPPERCASE with underscores
  - Examples: `VALID_MODALITIES`, `VALID_INTEGRATIONS`, `KNOWN_PLATFORMS`, `_DEFAULT_OBS_COLS`

## Where to Add New Code

**New Feature (analysis method):**
- Primary code: `sc_tools/tl/` → either directly in `sc_tools/tl/<method>.py` or in a sub-package like `sc_tools/tl/celltype/`
- Tests: `sc_tools/tests/test_<method>.py` (synthetic) and `test_<method>_real_data.py` (optional, with real data)
- Example: Adding RNA velocity analysis
  - Create: `sc_tools/tl/velocity.py` with function `run_rna_velocity(adata: AnnData, ...)`
  - Add to: `sc_tools/tl/__init__.py` → import and export
  - Test: `sc_tools/tests/test_velocity.py`

**New Preprocessing Step:**
- Primary code: `sc_tools/pp/<step>.py` (if large) or inline in `recipes.py` (if small)
- Integrate: Update `preprocess()` recipe in `sc_tools/pp/recipes.py` to call new step
- Tests: `sc_tools/tests/test_pp.py` (add fixture and test case)
- Example: Adding a new normalization method
  - Create: `sc_tools/pp/normalize.py` (exists; add function here)
  - Update: `sc_tools/pp/recipes.py` to conditionally call new normalization
  - Export: `sc_tools/pp/__init__.py`

**New Component/Module:**
- Implementation: `sc_tools/<component>/` directory with `__init__.py`, `*.py` files
- Public API: Define `__all__` list in `__init__.py`
- Parent export: Add import in `sc_tools/__init__.py` if module-level (e.g., `bm`, `memory`)
- Tests: `sc_tools/tests/test_<component>.py` (create file)
- Example: Adding a new benchmarking method
  - Create: `sc_tools/bm/new_method.py` with implementation
  - Export: Update `sc_tools/bm/__init__.py`
  - Test: Create `sc_tools/tests/test_new_benchmarking_method.py`

**Utilities:**
- Shared helpers: `sc_tools/utils/<category>.py`
- Private/internal helpers: prefix with `_` (e.g., `_local_path()` in loaders.py)
- Examples: `sc_tools/utils/signatures.py` (gene signature helpers), `sc_tools/utils/checkpoint.py` (save/load)

**HTML Reports/Templates:**
- Template files: `sc_tools/assets/<report_type>_template.html`
- Helper functions: `sc_tools/qc/report_utils.py` or `sc_tools/pl/benchmarking.py`
- Jinja2 rendering: Call `render_template()` from report_utils
- Example: Adding a new QC report stage
  - Create: `sc_tools/assets/<new_stage>_template.html` (Jinja2 template)
  - Function: `sc_tools/qc/report.py` → `generate_<new_stage>_report()`
  - Helper: `sc_tools/qc/report_utils.py` → add computation/rendering logic

**Database Migrations (optional):**
- Location: `sc_tools/migrations/versions/`
- Naming: `0NNN_<description>.py` (increment version number)
- Tool: Generated via `alembic revision --autogenerate -m "description"`
- Apply: `alembic upgrade head` during deployment
- Only needed if modifying registry ORM schema in `sc_tools/registry.py`

## Special Directories

**`sc_tools/assets/`:**
- Purpose: Static HTML/YAML templates for reports
- Generated: No (hand-written templates)
- Committed: Yes (version-controlled)
- Load method: `render_template(path, context)` in report_utils.py

**`sc_tools/data/`:**
- Purpose: Bundled biological resources (gene markers, benchmarks)
- Generated: Some files (IMC benchmark catalogs) are auto-generated; others are static
- Committed: Mostly yes; large benchmarks may be downloaded at runtime
- Load method: `sc_tools.data.io.load_marker_db()`, or direct file access

**`sc_tools/migrations/`:**
- Purpose: Database migration scripts (optional)
- Generated: Auto-generated by Alembic (`alembic revision --autogenerate`)
- Committed: Yes (required for reproducible schema upgrades)
- Only exists if SQLAlchemy/Alembic installed

**`sc_tools/tests/`:**
- Purpose: Unit and integration tests
- Generated: No (hand-written)
- Committed: Yes
- Run: `pytest sc_tools/tests/` or `pytest sc_tools/tests/test_pp.py` (specific)
- Real-data tests: `test_*_real_data.py` files are skipped unless real data fixture is available

---

*Structure analysis: 2026-03-20*
