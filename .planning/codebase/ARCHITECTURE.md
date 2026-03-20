# Architecture

**Analysis Date:** 2026-03-20

## Pattern Overview

**Overall:** Scanpy-inspired modular library with layered pipeline phases

**Key Characteristics:**
- Modality-agnostic (Visium, Visium HD, Xenium, CosMx, IMC compatible)
- Phase-based DAG orchestration (semantic slugs: `ingest_raw`, `qc_filter`, `preprocess`, `scoring`, `celltype_manual`, etc.)
- AnnData-centric data structure (in-memory with lazy storage backends)
- Registry-optional (database-free operation or optional SQLAlchemy/PostgreSQL for project tracking)
- GPU-aware with automatic fallback to CPU (rapids-singlecell vs scanpy)
- Pluggable storage backends (fsspec: local, S3, SFTP, GCS, Azure, Box)

## Layers

**Data Ingestion (Phase 0):**
- Purpose: Load raw vendor outputs (Space Ranger, Xenium Ranger, IMC output), validate manifests, generate HPC submission scripts
- Location: `sc_tools/ingest/`
- Contains: Platform-specific loaders (`load_visium_sample`, `load_xenium_sample`, `load_imc_sample`, etc.), batch manifest parsing, SLURM script generation
- Depends on: AnnData, scanpy, platform-specific reader libraries (spaceranger, tifffile, pillow)
- Used by: CLI scripts in `scripts/`, HPC-orchestrated pipelines

**Quality Control (Phases 1-2):**
- Purpose: Calculate QC metrics, filter low-quality cells, classify samples pass/fail, attach clinical metadata
- Location: `sc_tools/qc/`
- Contains: Metrics calculation (`calculate_qc_metrics`), spatial analysis (`spatially_variable_genes`), HTML report generation, per-sample QC classification
- Depends on: AnnData, statsmodels, squidpy, matplotlib
- Used by: Phase 1 QC filtering, Phase 2 metadata attachment, report generation

**Preprocessing & Integration (Phase 3):**
- Purpose: Normalize data, integrate batches, reduce dimensionality, cluster
- Location: `sc_tools/pp/`
- Contains: Normalization (log-transform, arcsinh for IMC), batch integration (scVI, Harmony, CytoVI, scanorama, BBKNN), dimensionality reduction (PCA, UMAP), clustering (Leiden), spatial clustering (UTAG)
- Depends on: scanpy, scipy, scikit-learn, optional integration libraries (scvi-tools, harmonypy, etc.)
- Used by: Phase 3 preprocessing recipe

**Analysis Tools (Phases 3.5+):**
- Purpose: Gene set scoring, statistical testing, colocalization, deconvolution, cell-type annotation
- Location: `sc_tools/tl/`
- Contains: Gene set management (`load_hallmark`, `load_msigdb_json`), enrichment testing (`run_ora`, `run_gsea_pseudobulk`), spatial colocalization, deconvolution (cell2location, tangram, destvi), cell-type annotation (sctype, celltypist, scArches, Geneformer, scGPT)
- Depends on: Optional model libraries (scvi-tools, tangram-sc, gseapy, decoupleR, celltypist, scarches, transformers)
- Used by: Phase 3.5+ analysis workflows

**Plotting & Visualization:**
- Purpose: Generate QC plots, spatial plots, heatmaps, benchmark comparison figures
- Location: `sc_tools/pl/`
- Contains: Wrappers around scanpy/matplotlib, spatial visualization, statistical annotations, benchmark plots (segmentation, integration), Plotly-based interactive reports
- Depends on: matplotlib, seaborn, scanpy, optional plotly
- Used by: Report generation, exploratory analysis

**Benchmarking:**
- Purpose: Compare segmentation and batch correction methods, generate metrics
- Location: `sc_tools/bm/`
- Contains: Segmentation scorers (detection, morphology, panoptic quality), integration metrics (scib-style), mask I/O adapters, SLURM orchestration, strategy runners (Cellpose, StarDist, DeepCell, HuggingFace models, SegFormer training)
- Depends on: scikit-image, scib-metrics, cellpose, stardist, deepcell, torch, torchvision, segmentation-models-pytorch
- Used by: Benchmarking pipelines, method comparison workflows

**Data Management:**
- Purpose: Storage abstraction, biological data platform registry, AnnData checkpoint validation
- Location: `sc_tools/storage.py`, `sc_tools/biodata.py`, `sc_tools/validate.py`
- Contains: fsspec URI resolution (`resolve_fs`, `open_file`, `read_h5ad`, `write_h5ad`), PlatformSpec registry (platform_for_project, get_platform, list_platforms), checkpoint validation (required obs/obsm per phase)
- Depends on: fsspec, pandas, AnnData
- Used by: All modules for I/O operations

**Pipeline Registry & Orchestration:**
- Purpose: Track project status, phases, checkpoints, provenance; provide MCP server for Claude integration
- Location: `sc_tools/registry.py`, `sc_tools/pipeline.py`, `sc_tools/mcp/`
- Contains: SQLAlchemy ORM (optional), phase DAG definition, checkpoint tracking, MCP tools server for phase status queries
- Depends on: SQLAlchemy (optional), alembic (optional), mcp (optional)
- Used by: HPC orchestrators, UI dashboards, phase-based automation

**Utilities:**
- Purpose: Logging, GPU detection, gene signature curation, figure saving
- Location: `sc_tools/utils/`, `sc_tools/memory/`
- Contains: Gene signature helpers, memory profiling, GPU environment detection, file saving utilities
- Depends on: Various (logging, psutil, torch for GPU detection)
- Used by: All modules

## Data Flow

**Phase 0a → Phase 0b (Ingestion):**

1. Raw vendor output directory (Space Ranger, Xenium, IMC, CosMx)
2. Load batch manifest (`load_batch_manifest`) - CSV or YAML with sample metadata
3. Validate manifest structure and paths (`validate_manifest`)
4. For each sample: call modality-specific loader (`load_visium_sample`, `load_xenium_sample`, etc.)
5. Each loader returns AnnData with standardized obs keys: `sample`, `library_id`, `raw_data_dir`, and obsm['spatial']
6. Concatenate samples: `concat_samples()` → single AnnData with obs['library_id'] batch key
7. Output: `results/adata.raw.h5ad` (Phase 0b checkpoint)

**Phase 1 (QC Filtering):**

1. Input: `adata.raw.h5ad` from Phase 0b
2. Calculate QC metrics: `calculate_qc_metrics()` adds obs columns (n_genes, n_counts, pct_mt, pct_hb)
3. Spatially variable genes: `spatially_variable_genes()` adds var columns (morans_i, etc.) if spatial
4. Classify samples: `compute_sample_metrics()` → per-sample pass/fail via `classify_samples()`
5. Filter cells: `filter_cells()` removes low-quality cells based on thresholds
6. Generate QC report: `generate_pre_filter_report()` with per-sample metrics tables and plots
7. Output: `results/adata.filtered.h5ad` (Phase 1 checkpoint)

**Phase 2 (Metadata Attachment):**

1. Input: `adata.filtered.h5ad` from Phase 1
2. Validate existing metadata structure
3. Attach clinical metadata: manual addition of obs columns (age, disease, treatment, etc.)
4. Regenerate post-filter QC report: `generate_post_filter_report()` pre vs post comparison
5. Output: `results/adata.metadata.h5ad` (Phase 2 checkpoint)

**Phase 3 (Preprocessing & Integration):**

1. Input: `adata.metadata.h5ad` from Phase 2
2. Call `pp.preprocess(adata, modality=..., batch_key="library_id", integration="scvi", ...)`
3. Orchestration steps:
   - Backup raw counts: `backup_raw()`
   - Normalize: `normalize_total()` or `normalize_imc()` (modality-aware)
   - Log transform: `log_transform()` or skip for IMC
   - Filter genes by pattern: `filter_genes_by_pattern()` (MT, RP, HB exclusions)
   - Highly variable gene selection (scanpy)
   - Batch integration (dispatch on `integration` param):
     - "scvi": `run_scvi()` → latent space in obsm['X_scvi']
     - "harmony": `run_harmony()` → corrected X + obsm['X_harmony']
     - "cytovi": `run_cytovi()` (IMC-only) → obsm['X_cytovi']
     - "none": skip integration
   - Dimensionality reduction: `pca()`, `neighbors()`, `umap()`
   - Clustering: `leiden()` or `cluster()` wrapper (resolution param)
   - Optional spatial clustering: if spatial_clustering="utag", `run_utag()`
4. Strategy selection (`select_strategy()`) for large vs small datasets affects computation resource allocation
5. GPU detection (`use_gpu="auto"`) enables rapids-singlecell if available
6. Generate post-integration QC report: `generate_post_integration_report()` with UMAP, cluster distribution
7. Output: `results/adata.integrated.h5ad` (Phase 3 checkpoint)

**Phase 3.5 (Scoring):**

1. Input: `adata.integrated.h5ad` from Phase 3
2. Gene set scoring: `score_signature()` for Hallmark, custom gene sets → obs columns with scores
3. Output: `results/adata.scored.h5ad` (Phase 3.5 checkpoint)

**Phase 4 (Cell Type Manual Annotation):**

1. Input: `adata.scored.h5ad` from Phase 3.5
2. Automated cell typing: `annotate_celltypes()` dispatcher across sctype, celltypist, scArches, Geneformer, scGPT, custom gates, ensemble
3. Manual curation: human expert adds obs['celltype_manual'] via UI or spreadsheet
4. Apply consensus: `apply_celltype_map()` merges automated predictions with manual input
5. Generate post-celltyping report: `generate_post_celltype_report()` with celltype abundance, UMAP by celltype
6. Output: `results/adata.celltypes.h5ad` (Phase 4 checkpoint)

**State Management:**

- **AnnData in-memory operations:** All transformations modify `adata.X`, `adata.obs`, `adata.obsm`, `adata.var`, `adata.uns` in-place unless `copy=True`
- **Checkpoint files:** Each phase outputs a standardized .h5ad checkpoint stored in `results/` with semantic naming
- **Registry tracking:** Optional SQLAlchemy backend records phase start/end times, file paths, validation status
- **Provenance:** uns['phase_history'] (or registry DB) logs which functions were called and with which parameters

## Key Abstractions

**AnnData-centric:** All functions accept/return AnnData objects. No custom data structures.

**Modality abstractions:**
- Purpose: Encapsulate platform-specific logic (normalization, integration, clustering parameters)
- Examples: `modality="visium"` vs `modality="imc"` (different normalization), `modality="visium_hd"` (single-cell resolution)
- Pattern: In `recipes.py` and loaders, `if modality == "imc": ...` branches for IMC-specific handling
- Location: `sc_tools/pp/recipes.py`, `sc_tools/ingest/loaders.py`, `sc_tools/qc/sample_qc.py`

**Integration strategies:**
- Purpose: Pluggable batch correction methods
- Examples: `run_scvi()`, `run_harmony()`, `run_cytovi()`, `run_combat()`, `run_bbknn()`, `run_scanorama()`, `run_scanvi()`, `run_resolvi()`
- Pattern: Each function returns the integration result in obsm['X_<method>'], optionally corrects X in-place
- Location: `sc_tools/pp/integrate.py`

**Cell type annotation methods:**
- Purpose: Pluggable automated annotation via different algorithms
- Examples: ScType (marker-based), CellTypist (ML-based), scArches (transfer learning), Geneformer (foundation model), scGPT (foundation model), custom gates (manual rules), ensemble (voting)
- Pattern: Each method returns obs['celltype_<method>'], dispatcher selects one or multiple
- Location: `sc_tools/tl/celltype/` with `_base.py`, `_sctype.py`, `_celltypist.py`, `_scarches.py`, `_geneformer.py`, `_scgpt.py`, `_custom_gates.py`, `_ensemble.py`

**Storage backends:**
- Purpose: Abstract filesystem API for local/cloud storage
- Examples: Local filesystem (`/path`), S3 (`s3://bucket/key`), SFTP (`sftp://host/path`), GCS (`gs://bucket/key`), Azure (`az://container/path`)
- Pattern: `resolve_fs(uri) → (AbstractFileSystem, path)`, then use fs.open(), fs.ls(), etc.
- Location: `sc_tools/storage.py`

**Pipeline phases:**
- Purpose: Define phase graph with dependencies and checkpoints
- Examples: PhaseSpec dataclass with label, depends_on, checkpoint, required_obs, required_obsm
- Pattern: DAG querying via `get_dag()`, `get_available_next(completed_phases)`, `get_phase_checkpoint(phase_slug)`
- Location: `sc_tools/pipeline.py`

## Entry Points

**CLI:**
- Location: `sc_tools/__main__.py`
- Triggers: `python -m sc_tools registry status`
- Responsibilities: Registry status display (database initialization, phase tracking query)

**Phase scripts (executed via HPC):**
- Location: `scripts/run_*phase*.py` (custom per-project, not in library)
- Triggers: Manual dispatch, HPC scheduler, orchestrator MCP tools
- Responsibilities: Load checkpoint from previous phase, call sc_tools functions, save checkpoint, update registry

**MCP server (Claude integration):**
- Location: `sc_tools/mcp/tools_server.py`, `sc_tools/mcp/registry_server.py`
- Triggers: Called by Claude Code orchestrator
- Responsibilities: Query registry status, update phase completion, provide phase DAG information

**Batch manifest parsing:**
- Location: `sc_tools/ingest/config.py`
- Triggers: `collect_all_batches()`, `load_batch_manifest()`
- Responsibilities: Parse CSV/YAML manifests, validate structure, return list of sample metadata dicts

## Error Handling

**Strategy:** Graceful degradation with detailed logging and optional validation

**Patterns:**

1. **Import-time guards:** Optional dependencies wrapped in try/except at module level
   - Example: `sc_tools/registry.py` imports SQLAlchemy only if installed
   - Example: `sc_tools/pp/_gpu.py` detects RAPIDS availability and falls back to CPU scanpy

2. **Validation functions:** Explicit checkpoint validation via `validate.py`
   - `validate_checkpoint(adata, phase="qc_filter")` returns list of validation issues (not raising)
   - `validate_file(filepath, phase="p1")` reads file and validates

3. **Warnings for deprecated features:**
   - Example: Legacy phase codes (p1, p2, p3, p35, p4) emit DeprecationWarning and map to semantic slugs
   - Example: Old-style checkpoint filenames (adata.raw.h5ad, adata.filtered.h5ad) still accepted

4. **Logging at INFO/DEBUG levels:**
   - Example: `logger.info("Loaded Visium sample %s: %d spots x %d genes", sample_id, n_obs, n_vars)`
   - Example: `logger.debug("Starting batch integration with method: %s", integration_method)`

5. **Explicit error raising:**
   - Unknown modality: `raise ValueError(f"Unknown modality: {modality}. Choose from {VALID_MODALITIES}")`
   - Missing required obs: `raise ValueError(f"Missing required obs column: {col}")`
   - Storage backend not installed: `raise ImportError("fsspec backend for 's3://' not available. Install: pip install s3fs")`

## Cross-Cutting Concerns

**Logging:**
- Framework: Python logging via `logger = logging.getLogger(__name__)`
- Configuration: pyproject.toml specifies test filterwarnings (ignore DeprecationWarning, FutureWarning)
- Pattern: Info level for major steps (sample load, integration start), Debug level for details

**Validation:**
- Checkpoint validation via `validate.py` (required obs/obsm per phase)
- Manifest validation via `ingest/config.py` (required columns, file existence)
- Function-level parameter validation (modality in VALID_MODALITIES, etc.)

**Authentication:**
- Cloud storage (S3, GCS, Azure) credentials: handled by fsspec/environment variables
- Box storage: OAuth flow at runtime (optional dependency)
- Registry DB: SQLAlchemy connection string (optional, defaults to SQLite)

**Caching:**
- AnnData `.h5ad` files serve as checkpoints (persistent cache)
- Optional in-memory caching in integration algorithms (scVI, Harmony, etc.)
- No explicit application-level caching layer

---

*Architecture analysis: 2026-03-20*
