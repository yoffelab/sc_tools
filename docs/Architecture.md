# Project Architecture: Scalable Spatial Omics Analysis

This document outlines the directory structure, data flow, and script inventory. The layout is **scalable across projects and data types** (Visium, Visium HD, Xenium, IMC, CosMx). All **metadata**, **results**, and **figures** are **project-specific** under `projects/<platform>/<project_name>/`. There is no repo-root `metadata/`, `results/`, or `figures/`. The pipeline is non-linear; see `README.md` (Pipeline Workflow section) for the diagram.

## 1. Directory Overview

```text
.
├── CLAUDE.md               # Claude Code project config (repo root)
├── README.md               # GitHub landing page (repo root)
├── pyproject.toml          # Package build (sc_tools installable)
├── environment.yml         # Conda environment
├── requirements.txt        # pip dependencies
├── Makefile                # Pipeline orchestration (project-aware; default PROJECT=projects/visium/ggo_visium)
│
├── docs/                   # All documentation
│   ├── Architecture.md     # System roadmap (this file)
│   ├── Mission.md          # Toolkit todo list and roadmap; project-specific in projects/<type>/<name>/Mission.md
│   ├── Journal.md          # Repo-level decision log; project-specific in projects/<type>/<name>/Journal.md
│   ├── journal_summary.md  # Short summary of Journal.md; per-project under projects/<type>/<name>/
│   ├── skills.md           # Mandatory coding and statistical standards
│   ├── project_setup.md    # Environment and container setup notes
│   ├── Scratch Pad.md      # Exploratory notes
│   ├── wiki/               # Obsidian vault (symlinks + wiki-native + .gen.md files)
│   └── ...                 # Sphinx source (conf.py, index.rst, api/, tutorials/)
│
├── sc_tools/               # Reusable Python package (scanpy-style API) — NOT project-specific
│   ├── pl/                  # Plotting: spatial, heatmaps, statistical, volcano, save
│   ├── tl/                  # Tools: testing, colocalization, deconvolution, io, celltype (automated cell typing)
│   ├── qc/                  # QC: calculate_qc_metrics, filter_cells, filter_genes, highly_variable_genes, spatially_variable_genes
│   ├── pp/                  # Preprocessing: normalization, integration, clustering, recipes
│   ├── ingest/              # ingest_raw: batch manifests, command builders, modality loaders
│   ├── bm/                  # Benchmarking: integration scoring, IMC segmentation benchmark, CLI, reports
│   ├── data/                # Bundled reference data: Hallmark signatures, IMC benchmark data, QC/report templates
│   ├── pipeline.py          # Phase DAG: PhaseSpec, STANDARD_PHASES, get_available_next(), get_phase_checkpoint()
│   ├── storage.py           # fsspec URI resolution + smart read/write (local, S3, SFTP, GCS, Box)
│   ├── registry.py          # SQLAlchemy registry: projects, datasets, SLURM jobs, agent tasks, project_phases
│   ├── migrations/          # Alembic migration scripts for registry schema evolution
│   ├── mcp/                 # FastMCP servers: sc-tools (analysis tools) + sc-registry (bookkeeping)
│   ├── validate.py          # Checkpoint validation (p1-p4) per Architecture.md Section 2.2
│   ├── memory/              # Profiling, GPU
│   ├── utils/               # Signatures, versioned save helpers
│   └── tests/               # Package unit tests (pytest)
│
├── scripts/                # Shared/legacy scripts; project scripts live under projects/<platform>/<project>/scripts/
│
├── projects/               # All project-specific content lives here
│   ├── create_project.sh   # Create new project: ./projects/create_project.sh <project_name> <data_type>
│   ├── README.md           # How to create projects and use Mission/Journal
│   ├── visium/
│   │   └── ggo_visium/     # Example project
│   │       ├── data/       # Raw sequencing and imaging
│   │       ├── figures/    # QC/raw, QC/post, manuscript/, etc.
│   │       ├── metadata/   # gene_signatures.json, sample_metadata.csv, celltype_map.json, obs.csv
│   │       ├── scripts/    # Project-specific analysis scripts
│   │       ├── results/    # AnnData (.h5ad), CSVs
│   │       ├── outputs/    # Intermediate outputs (deconvolution logs, etc.)
│   │       ├── tests/      # Project integration tests (pytest)
│   │       ├── Mission.md
│   │       ├── Journal.md
│   │       └── journal_summary.md
│   ├── visium_hd/
│   ├── xenium/             # Placeholder — no active project yet (.gitkeep only)
│   ├── imc/
│   ├── cosmx_1k/
│   ├── cosmx_6k/
│   └── cosmx_full_library/
└── .gitignore, .githooks/
```

**Creating a new project:** Run `./projects/create_project.sh <project_name> <data_type>`. Valid `data_type`: `visium` | `visium_hd` | `visium_hd_cell` | `xenium` | `imc` | `cosmx_1k` | `cosmx_6k` | `cosmx_full_library`.

---

## 2. Project-Specific Paths and Checkpoint Nomenclature

All of the following live under `projects/<platform>/<project_name>/`. **Checkpoint AnnData files must use the standard names below** so tooling can validate metadata and downstream scripts can assume consistent paths.

**Phase DAG:** Pipeline phases are defined as semantic slugs in `sc_tools/pipeline.py` (`STANDARD_PHASES`). Each slug has explicit `depends_on` references so tooling can discover available next steps and validate dependency order. Use `get_available_next(completed)` to query what can run next, and `get_phase_checkpoint(slug)` to retrieve the expected output path. New branches or project-specific phases can be registered via `extend_dag(slug, PhaseSpec(...))` without changing library code.

**Phase slug → old code mapping:**

<!-- PHASE_TABLE:START -->

| Slug | Old code | Name | Checkpoint | Required Data | QC Report |
|------|----------|------|------------|---------------|-----------|
| `ingest_raw` | p0a | Raw Data Processing | - | - |  |
| `ingest_load` | p0b | Load into AnnData | `data/{sample_id}/adata.ingested.h5ad` | `obs[sample, library_id, raw_data_dir]`, `obsm[spatial]`, `X` raw counts |  |
| `qc_filter` | p1 | QC Filtering + Concatenation | `results/adata.filtered.h5ad` | `obs[sample, raw_data_dir]`, `obsm[spatial]`, `X` raw counts, concatenated | `pre_filter_qc_{date}.html` |
| `metadata_attach` | p2 | Metadata Attachment | `results/adata.annotated.h5ad` | `obs[sample, raw_data_dir]`, `obsm[spatial]`, `X` raw counts, concatenated | `post_filter_qc_{date}.html` |
| `preprocess` | p3 | Normalize + Integrate + Cluster | `results/adata.normalized.h5ad` | `obs[leiden]`, `obsm[X_scvi]`, `X` normalized (adata.raw backed up) | `post_integration_qc_{date}.html` |
| `demographics` | p3.5 | Cohort Demographics | - | - |  |
| `scoring` | p3.5b | Gene Scoring + Auto Cell Typing | `results/adata.scored.h5ad` | `obsm[signature_score, signature_score_z]`, `X` normalized |  |
| `celltype_manual` | p4 | Manual Cell Typing | `results/adata.celltyped.h5ad` | `obs[celltype, celltype_broad]`, `X` normalized | `post_celltyping_qc_{date}.html` |
| `biology` | p5 | Downstream Biology | - | - |  |
| `meta_analysis` | p6/p7 | Meta Analysis | - | - |  |

<!-- PHASE_TABLE:END -->

### 2.1 Mandatory checkpoint filenames (results/)

| Slug | Standard path | Description |
|------|----------------|-------------|
| `ingest_raw` | `data/{sample_id}/outs/` | Per-sample raw platform output (Space Ranger, Xenium Ranger, IMC). |
| `ingest_load` | `data/{sample_id}/adata.ingested.h5ad` | Per-sample AnnData loaded from platform output via `sc_tools.ingest.loaders`. |
| `ingest_load` | `data/{sample_id}/spatialdata.zarr` | Per-sample SpatialData object (optional; Visium HD / Xenium when rich spatial data needed). |
| `ingest_raw` | `metadata/phase0/all_samples.tsv` | Collected batch manifest (auto-generated by `collect_all_batches()`). |
| `qc_filter` | `results/adata.filtered.h5ad` | Concatenated, QC-filtered AnnData across all samples. |
| `metadata_attach` | `results/adata.annotated.h5ad` | AnnData with clinical/sample metadata (and masks if applicable) attached. |
| `preprocess` | `results/adata.normalized.h5ad` | Filtered, normalized, batch-corrected, clustered (no cell typing). |
| `scoring` | `results/adata.scored.h5ad` | Normalized + gene signature scores (and automated cell typing if applied). |
| `celltype_manual` | `results/adata.celltyped.h5ad` | Manual cell type labels applied (celltype, celltype_broad). |
| `biology`/`meta_analysis` (aggregated) | `results/adata.{level}.{feature}.h5ad` | Aggregated by `level` in `roi` or `patient` and `feature` in `mean_expression` or `celltype_frequency`. |

**Legacy names (migration):** Existing projects may still use `adata.p0.h5ad`, `adata.raw.h5ad`, `adata.raw.p1.h5ad`, `adata.annotated.p2.h5ad`, `adata.normalized.p3.h5ad`, `adata.normalized.scored.p35.h5ad`, `adata.celltyped.p4.h5ad`. Both old and new names are accepted during transition; `validate_checkpoint()` reports which convention is in use.

**Aggregated examples:** `results/adata.roi.mean_expression.h5ad`, `results/adata.patient.celltype_frequency.h5ad`.

### 2.2 Required metadata (for validation)

Checkpoints must satisfy the following so scripts and validators can rely on them. Tooling (e.g. `sc_tools` or a project test) should **check** these before using a checkpoint.

| Checkpoint | Required `adata` contents |
|------------|---------------------------|
| **adata.ingested.h5ad** (per-sample, `ingest_load`) | `obs['sample']`, `obs['library_id']`, `obs['raw_data_dir']`, `obsm['spatial']`; `X` raw counts; single sample only (not concatenated). Set by `sc_tools.ingest.loaders`. **IMC only (optional):** `adata.uns['spatial'][sample_id]` may contain image data when `load_imc_sample(load_images=True)` is used — see IMC image schema below. |
| **adata.filtered.h5ad** (`qc_filter`) | `obs['sample']`, `obs['raw_data_dir']` (or equivalent), `obsm['spatial']`; `X` raw counts; no normalization; all samples concatenated. |
| **adata.annotated.h5ad** (`metadata_attach`) | All of `qc_filter`; `obs` includes clinical/metadata columns from `metadata/sample_metadata.csv` (or join equivalent); optional `uns` keys for images/masks. |
| **adata.normalized.h5ad** (`preprocess`) | Normalized/batch-corrected representation (e.g. `obsm['X_scvi']`); `obs['leiden']` (or cluster column); `adata.raw` backed up. |
| **adata.scored.h5ad** (`scoring`) | All of `preprocess` (or `metadata_attach` where preprocess not used); signature scores in `obsm['signature_score']` (raw) and `obsm['signature_score_z']` (z-scored), column names = full path (e.g. `Myeloid/Macrophage_Core`); `uns['signature_score_report']` for per-signature n_present, n_missing, status; optional `obs['celltype']`/`celltype_broad` if automated typing run. |
| **adata.celltyped.h5ad** (`celltype_manual`) | All of `scoring`; `obs['celltype']`, `obs['celltype_broad']` from `metadata/celltype_map.json`. |
| **adata.{level}.{feature}.h5ad** (`meta_analysis`) | `obs` indexed by `level` (roi or patient); `X` or layer holds aggregated `feature` (mean expression or celltype frequency). |

### 2.2b IMC image schema in adata.uns['spatial'] (optional `ingest_load` phase)

When `load_imc_sample(load_images=True)` is called, the per-ROI TIFF stack is loaded and stored in the standard squidpy/scanpy format, making `sc.pl.spatial(img_key="hires")` work without modification:

```
adata.uns['spatial'][sample_id]
  images/
    hires   (H, W, 3) uint8    — percentile-clipped RGB composite (PanCK=R, CD3=G, DNA1=B default)
    full    (C, H, W) float32  — arcsinh(x/5) normalized full channel stack
    mask    (H, W) int32       — labeled cell segmentation mask (if load_mask=True)
  scalefactors/
    tissue_hires_scalef: 1.0/downsample
    tissue_lowres_scalef: 1.0/downsample
    spot_diameter_fullres: 1.0   (1 pixel ~ 1 um in IMC)
  metadata/
    channels: list[str]         — ordered protein names (before "(")
    channel_strings: list[str]  — full MarkerName(IsotopeTag) strings
    rgb_channels: {R, G, B}     — resolved protein names used for composite
    rgb_indices: {R, G, B}      — TIFF stack channel indices used
    pixel_size_um: 1.0
```

**TIFF file naming** (ElementoLab/steinbock convention):
- `processed/{sample}/tiffs/{roi_id}_full.tiff` — (C, H, W) multi-channel intensity stack
- `processed/{sample}/tiffs/{roi_id}_full.csv` — channel index CSV (`index, MarkerName(IsotopeTag)`)
- `processed/{sample}/tiffs/{roi_id}_full_mask.tiff` — labeled cell segmentation mask
- `processed/{sample}/tiffs/{roi_id}_full_nucmask.tiff` — labeled nuclear mask
- `processed/{sample}/tiffs/{roi_id}_Probabilities.tiff` — ilastik probability map

**Channel CSV format** (`*_full.csv`): two columns, first = 0-based index, second = `MarkerName(IsotopeTag)` (e.g. `CD3(Er170)`, `PanCK(Pt195)`, `<EMPTY>(In115)`).

**Panel CSV format** (`channel_labels.csv` / `data/{sample}.channel_labels.csv`): columns `channel, Target, Metal_Tag, Atom, full, ilastik`. `Target` = protein name; `ilastik=1` = used for ilastik segmentation.

**`IMCPanelMapper`** (`sc_tools.ingest.IMCPanelMapper`): resolves any name form (protein name, `MarkerName(IsotopeTag)`, isotope tag, partial substring) to a TIFF stack channel index. Reads `*_full.csv` via `from_full_csv()` and optionally `channel_labels.csv` via `from_panel_csv()`.

### 2.3 Other project paths (metadata, optional results, figures)

| Path | Description |
|------|-------------|
| `metadata/sample_metadata.csv` or `.xlsx` | Sample→clinical metadata map. Enables `metadata_attach` bypass. |
| `metadata/celltype_map.json` | cluster_id→celltype mapping for `celltype_manual`. |
| `metadata/gene_signatures.json` | Gene signatures for scoring (and per-signature `metadata/{name}.json` for obsm storage). |
| `results/adata.deconvolution.h5ad` | Cell-type proportions (optional; `scoring`). |
| `results/tmp/integration_test/{method}.h5ad` | Intermediate integration test results from benchmark subsample. Kept for re-annotation and further scoring without re-running integration. |
| `results/integration_method.txt` | Records which integration method was selected and used for the final dataset. |

**figures/ subdirectory layout** — all projects use this 6-directory structure (created by `create_project.sh`). The hook at `~/.claude/hooks/figure_quality_check.py` enforces different evaluation standards per subdirectory.

| Directory | Hook standard | Purpose |
|-----------|---------------|---------|
| `figures/scratch/` | **Skip** — never evaluated | Throwaway exploratory plots; use freely without triggering review |
| `figures/QC/` | **Technical only** (Part 1) | Auto-generated pipeline QC reports (`pre_filter_qc_YYYYMMDD.html`, `post_filter_qc_YYYYMMDD.html`, `post_integration_qc_YYYYMMDD.html`, `post_celltyping_qc_YYYYMMDD.html`) |
| `figures/exploratory/` | **Standard** (Parts 1 + 2) | Analysis in progress — real figures not yet curated for manuscript |
| `figures/deconvolution/` | **Standard** (Parts 1 + 2) | Cell-type deconvolution outputs (`scoring`); per-method subdirs (`cell2location/`, `tangram/`) are recognised automatically |
| `figures/insightful/` | **High** (Parts 1 + 2, high bar) | Clear biological findings; candidate figures for manuscript or key communication |
| `figures/manuscript/` | **Strict** (Parts 1 + 2, publication bar) | Final publication figures — must be fully publication-ready |
| `figures/supplementary/` | **Strict** (Parts 1 + 2, publication bar) | Supplementary material — same standard as manuscript |

Files named `*_draft.*`, `*_test.*`, `*_tmp.*`, or `*_temp.*` are skipped by the hook regardless of directory. The hook writes its evaluation to `{figure_stem}_claude_analytics_and_legend.txt` alongside each figure, starting with `SATISFACTORY: YES|NO` and `REVISION: N`. Re-evaluation is triggered automatically until `SATISFACTORY: YES` (max 5 attempts before human escalation).

**Legacy / migration:** Existing projects may still use `adata.annotation.masked.h5ad`, `scvi.leiden.phenotyped.h5ad`, `adata.img.genescores.h5ad` until scripts are updated. New pipelines and scripts **must** write the standard checkpoint names above. Validators should accept either legacy or standard names and report which convention is used.

**Snakemake:** Use `$(PROJECT)` or config-based paths for all project paths, e.g. `$(PROJECT)/results/adata.filtered.h5ad`.

### 2.4 Signature scoring (`scoring`): how to run, inputs, outputs

Gene signature scoring uses `sc_tools.tl.score_signature` and writes **obsm** (`signature_score`, `signature_score_z`) and **uns** (`signature_score_report`). Each project has a thin script; run via Snakemake from **project root** (so paths like `results/` and `metadata/` resolve).

| Project | How to run | Inputs | Output |
|--------|------------|--------|--------|
| **Robin** | `snakemake -d projects/visium_hd/robin -s projects/visium_hd/robin/Snakefile scoring` or `cd projects/visium_hd/robin && python scripts/run_signature_scoring.py` | `results/adata.normalized.h5ad`, `metadata/gene_signatures.json` | `results/adata.scored.h5ad` |
| **ggo_visium** | `snakemake -d projects/visium/ggo_visium -s projects/visium/ggo_visium/Snakefile scoring` or `cd projects/visium/ggo_visium && python scripts/score_gene_signatures.py` | `results/adata.annotated.h5ad`, `metadata/gene_signatures.json` | `results/adata.scored.h5ad` |

**Note:** Robin builds scoring from **preprocess** (normalized adata). ggo_visium builds scoring from **metadata_attach** (annotated adata). Downstream scripts in both projects read the scoring checkpoint.

**Snakemake + containers:** The pipeline uses Snakemake as the workflow engine with Apptainer/Singularity (Linux/HPC, primary) and Docker (macOS/Windows, fallback). Auto-configure via `scripts/run_container.sh`; publish container image(s) to a registry for HPC pull. See Mission.md CI/CD Roadmap.

### 2.5 QC Report Workflow (4 date-versioned HTML reports)

The pipeline produces four QC HTML reports at key checkpoints. Each report is date-versioned (`YYYYMMDD`) for auditability.

| # | Report | Phase | Primary Metrics | Notes |
|---|--------|-------|-----------------|-------|
| 1 | `pre_filter_qc_YYYYMMDD.html` | `qc_filter` (entry) | Per-sample cell/spot counts, gene/protein detection, %MT | Raw data overview before any filtering. |
| 2 | `post_filter_qc_YYYYMMDD.html` | `qc_filter`/`metadata_attach` (exit) | Pre-vs-post comparison, HVG/SVG, sample pass/fail | After QC filtering and metadata attachment. |
| 3 | `post_integration_qc_YYYYMMDD.html` | `preprocess` (exit) | **Batch score (primary)**, UMAP grid, cluster distribution, integration benchmark | After integration. Bio metrics (ARI, NMI, ASW celltype) are informational only — celltype labels are preliminary (leiden clustering). |
| 4 | `post_celltyping_qc_YYYYMMDD.html` | `celltype_manual` (exit) | **Full bio + batch scores**, re-scored integration methods | After validated celltyping. Bio conservation metrics are now meaningful. Re-evaluates all candidate integrations if `results/tmp/integration_test/` exists. |

**Key principle:** Batch score is the primary metric for integration method selection (report 3). Bio metrics become the primary evaluation criterion only after celltyping is validated (report 4). This avoids circular reasoning where integration is optimized for celltype labels derived from that same integration.

---

## 3. Phase Details and Data Flow

### ingest_raw: Upstream Raw Data Processing

**Purpose:** Run platform-specific upstream tools on HPC and produce per-sample AnnData / SpatialData objects ready for QC and concatenation in `qc_filter`. `ingest_raw` has two sub-steps:

#### ingest_raw: Platform tools (HPC)

Run Space Ranger, Xenium Ranger, or the IMC pipeline on raw FASTQs / images to produce platform output directories.

**Batch manifest system:**
- Per-batch TSV files under `metadata/phase0/` (e.g., `batch1_samples.tsv`, `batch2_samples.tsv`)
- `sc_tools.ingest.collect_all_batches()` concatenates all batch TSVs into `metadata/phase0/all_samples.tsv`
- Each TSV has modality-specific required columns (see `sc_tools.ingest.config.REQUIRED_COLUMNS`)

| Modality | TSV columns |
|----------|-------------|
| Visium | `sample_id`, `fastq_dir`, `image`, `slide`, `area`, `batch` |
| Visium HD / Visium HD Cell | `sample_id`, `fastq_dir`, `cytaimage`, `slide`, `area`, `batch` |
| Xenium | `sample_id`, `xenium_dir`, `batch` |
| IMC | `sample_id`, `mcd_file`, `panel_csv`, `batch` |
| CosMx | `sample_id`, `cosmx_dir`, `panel_tier` (1k/6k/full_library), `batch` — flat files / RDS converted to AnnData in `ingest_load` |

**Snakemake:** `spaceranger_count` (or equivalent) per-sample rule; `phase0` target expands all samples from `all_samples.tsv`.

**Command builders:** `sc_tools.ingest.spaceranger.build_spaceranger_count_cmd()`, `sc_tools.ingest.xenium.build_xenium_ranger_cmd()`, `sc_tools.ingest.imc.build_imc_pipeline_cmd()`.

**Output:** `data/{sample_id}/outs/` per sample.

#### ingest_load: Load into AnnData / SpatialData

Convert per-sample platform outputs into portable AnnData (or SpatialData) objects with standardized metadata. This is the `ingest_raw` checkpoint consumed by `qc_filter`.

**Loaders (`sc_tools.ingest.loaders`):**

| Modality | Loader | Output format |
|----------|--------|---------------|
| Visium | `load_visium_sample()` | `data/{sample_id}/adata.ingested.h5ad` |
| Visium HD | `load_visium_hd_sample()` | `data/{sample_id}/adata.ingested.h5ad` or `data/{sample_id}/spatialdata.zarr` |
| Visium HD Cell | `load_visium_hd_cell_sample()` | `data/{sample_id}/adata.ingested.h5ad` |
| Xenium | `load_xenium_sample()` | `data/{sample_id}/adata.ingested.h5ad` or `data/{sample_id}/spatialdata.zarr` |
| IMC | `load_imc_sample()` | `data/{sample_id}/adata.ingested.h5ad` |
| CosMx | `load_cosmx_sample()` | `data/{sample_id}/adata.ingested.h5ad` |

**CosMx loading (deprioritized):** CosMx output is flat CSV / Parquet files (or RDS from NanoString software). `load_cosmx_sample()` reads the expression matrix, FOV coordinates, and cell metadata; converts to AnnData with `obsm['spatial']` set from cell centroid coordinates (x, y in microns). RDS inputs require `rpy2` + `anndata2ri`; flat file inputs use `pandas` + `squidpy` for spatial graph construction. Three panel tiers exist with different scale characteristics:

- **CosMx 1k:** ~1,000-plex targeted panel. Standard flat CSV/Parquet export from AtoMx.
- **CosMx 6k:** ~6,000-plex panel. Larger expression matrices; same flat file format.
- **CosMx full_library:** Whole-transcriptome (~18k genes). Significantly larger files; may require chunked loading or backed AnnData (`adata.X` as sparse on disk).

Each loader sets `obs['sample']`, `obs['library_id']`, `obs['raw_data_dir']`, `obsm['spatial']`, and ensures `X` contains raw counts. See Section 2.2 for full requirements.

**SpatialData (optional):** For Visium HD and Xenium, `spatialdata.zarr` is preferred when downstream analysis requires full image pyramids, masks, or subcellular coordinates. Use `spatialdata-io` for loading; the `.zarr` store can coexist with `adata.ingested.h5ad` for tools that require AnnData.

**Output:** `data/{sample_id}/adata.ingested.h5ad` (required) and/or `data/{sample_id}/spatialdata.zarr` (optional).

### Checkpoint Validation

Each checkpoint (p1 through p4) has a validation sentinel rule in the Snakefile. The sentinel file `results/.adata.{name}.validated` is created by `scripts/validate_checkpoint.py` which calls `sc_tools.validate.validate_checkpoint()`. Downstream rules depend on the sentinel, ensuring metadata contracts (Section 2.2) are enforced before the next phase.

**CLI:** `python scripts/validate_checkpoint.py <path> --phase <qc_filter|metadata_attach|preprocess|scoring|celltype_manual> [--fix] [--warn-only]`

Legacy phase codes (`p1`, `p2`, `p3`, `p35`, `p4`) are still accepted but emit a deprecation warning.

**Auto-fix (--fix):** Renames `obs['batch']` to `obs['raw_data_dir']` (p1 only). All other issues require human intervention.

### qc_filter: QC and Concatenation

**Input:** Per-sample `data/{sample_id}/adata.ingested.h5ad` objects produced in `ingest_load`.

**Steps:**
1. Load all per-sample `ingest_load` AnnData objects.
2. Apply per-sample QC: `sc_tools.qc.filter_spots()` (modality-aware thresholds); compute metrics via `sc_tools.qc.calculate_qc_metrics()`.
3. Classify and optionally remove low-quality samples: `sc_tools.qc.sample_qc.classify_samples()`, `apply_qc_filter()`.
4. Concatenate across samples: `sc_tools.ingest.concat_samples()`.
5. Generate QC reports: 2x2 histograms, violin metrics, spatial %MT plots, cross-sample comparison bar/violin/scatter. Pre-filter → `figures/QC/raw/`.

**Required annotations on output:** `obs['sample']`, `obs['raw_data_dir']`, `obsm['spatial']`; `X` raw counts; no normalization. See Section 2.2.

**QC (sc_tools.qc):** `calculate_qc_metrics`, `filter_cells`, `filter_genes`, `highly_variable_genes`; squidpy spatially variable genes. HTML report via `sc_tools.qc.report.generate_qc_report()`.

**Outputs:** `$(PROJECT)/results/adata.filtered.h5ad`, `$(PROJECT)/figures/QC/raw/*`. Must satisfy required metadata for `qc_filter` (Section 2.2).

---

### metadata_attach: Metadata Attachment (Human-in-Loop)

**Bypass:** Provide `$(PROJECT)/metadata/sample_metadata.csv` or `.xlsx`.

**Without file:** Human prepares map; cannot skip automatically.

**Outputs:** `$(PROJECT)/results/adata.annotated.h5ad`. Must satisfy required metadata for `metadata_attach` (Section 2.2).

---

### preprocess: Preprocessing

Backup `adata.raw`; filter; normalize; batch correct; cluster. Post-QC → `$(PROJECT)/figures/QC/post/`. No automated cell typing (that is in `scoring`).

#### Integration Benchmark Workflow (preprocess sub-step)

Before committing to a single integration method for the full dataset, run a benchmark on a subsample:

1. **Subsample:** Select ~20 samples (or all if fewer) to keep benchmark tractable.
2. **Run all candidate methods:** Each method produces an embedding stored in a temporary AnnData under `results/tmp/integration_test/{method}.h5ad`. Methods tested depend on modality:
   - **IMC:** Harmony, ComBat, scVI (raw), scANVI (raw), CytoVI (arcsinh), IMC Phenotyping, IMC Pheno+Harmony, Z-score+Harmony, Unintegrated PCA.
   - **Transcriptomic (Visium, Xenium, CosMx):** Harmony, ComBat, scVI, scANVI (if celltypes), resolVI (spatial-aware), BBKNN, Scanorama, Unintegrated PCA.
   - **Visium HD (bin-level, `visium_hd`):** Same as Transcriptomic.
   - **Visium HD Cell (`visium_hd_cell`):** Same as Transcriptomic. No deconvolution step — cell segmentation (SpaceRanger 4) provides single-cell resolution directly; `adata.deconvolution.h5ad` is not produced.
3. **Benchmark (batch-focused):** `compare_integrations()` computes batch metrics (ASW batch, PCR, graph connectivity) and informational bio metrics (ASW celltype, ARI, NMI — based on preliminary leiden clustering, not validated celltypes). **Batch score is the primary selection criterion** at this stage because celltype labels are not yet validated. Bio metrics are reported but should not drive method selection pre-celltyping.
4. **Select best method:** Automatically pick the method with the highest batch score. Record the choice in `results/integration_method.txt`.
5. **Integrate full dataset:** Apply the selected method to all samples. Save `results/adata.normalized.h5ad`.
6. **Generate post-integration QC report:** `figures/QC/post_integration_qc_YYYYMMDD.html` with batch score highlighted as the primary metric.

The temporary `results/tmp/integration_test/{method}.h5ad` files are kept for later use: after celltyping (`celltype_manual`), the same embeddings can be re-scored with validated celltypes to produce the post-celltyping QC report.

**Outputs:** `$(PROJECT)/results/adata.normalized.h5ad`. Must satisfy required metadata for `preprocess` (Section 2.2).

---

### demographics: Demographics (Branching)

Separate branch from `preprocess` (parallel to `scoring`). sc_tools helpers: piechart, histogram, violinplot, barplot, stacked barplot, scatterplot, correlogram, heatmap. Figure 1 for cohort description.

---

### scoring: Gene Scoring, Automated Cell Typing, Deconvolution

Separate branch from `preprocess` (parallel to `demographics`); connects to `celltype_manual`. Always apply basic gene sets (e.g. Hallmark, loaded via `sc_tools.tl.load_hallmark()`) and any project-provided signatures from `metadata/{signature_name}.json`. Store scores in **`adata.obsm['signature_score']`** (raw) and **`adata.obsm['signature_score_z']`** (z-scored) via `sc_tools.tl.score_signature`; column names are full paths (e.g. `Hallmark/HYPOXIA`, `Myeloid/Macrophage_Core`). Scores are **not** stored in `obs` by default. Automated cell typing (cluster → celltype). For Visium and Visium HD only: optional cell-type deconvolution (DestVI, Cell2location, Tangram) → `$(PROJECT)/results/adata.deconvolution.h5ad`.

**Deconvolution applicability by modality:**
- **Visium (`visium`):** Deconvolution applies — each spot covers ~10–50 cells.
- **Visium HD (`visium_hd`, bin-level):** Deconvolution applies — 8 µm bins still capture multiple cells.
- **Visium HD Cell (`visium_hd_cell`), Xenium, CosMx, IMC:** Deconvolution is **skipped** — all produce single-cell resolution data (cell segmentation or targeted single-cell capture). Cell type composition is determined directly from clustering and cell typing.

**Outputs:** `$(PROJECT)/results/adata.scored.h5ad` (must satisfy Section 2.2); optional `adata.deconvolution.h5ad`. Required for `biology`.

---

### celltype_manual: Manual Cell Typing (Human-in-Loop; Skippable)

Skippable if automated cell typing in `scoring` is adequate. JSON format `{cluster_id: {celltype_name, total_obs_count}}`; match cluster_id type; produce celltype and celltype_broad. Iterative until satisfactory. Save to `$(PROJECT)/metadata/celltype_map.json`.

#### Post-Celltyping QC Report (celltype_manual exit)

After celltyping is finalized, generate `figures/QC/post_celltyping_qc_YYYYMMDD.html`. This report re-evaluates integration quality using **validated celltype labels**, making bio conservation metrics (ARI, NMI, ASW celltype, isolated labels, cLISI) fully meaningful. The report includes:

- Full integration benchmark table with both batch and bio scores (bio metrics now trustworthy).
- Comparison to the post-integration report (which used preliminary leiden labels).
- Note on which integration method was used (from `results/integration_method.txt`).
- If `results/tmp/integration_test/{method}.h5ad` files exist, re-score all candidate methods with the validated celltypes to show whether the batch-score-based selection at `preprocess` also performs well on bio metrics. This reveals if a different method would have been better given the final celltypes.

**Outputs:** `$(PROJECT)/results/adata.celltyped.h5ad`. Must satisfy required metadata for `celltype_manual` (Section 2.2).

---

### biology: Downstream Biology

Uses gene scores and (optionally) deconvolution from `scoring`. Spatial/process analysis, colocalization, neighborhood enrichment, publication figures. Reads from `adata.scored.h5ad` (or `adata.celltyped.h5ad` if `celltype_manual` was run).

---

### meta_analysis: Meta Analysis (Optional)

Aggregate ROI/patient; downstream on aggregated data.

**Outputs:** `$(PROJECT)/results/adata.{level}.{feature}.h5ad` with `level` ∈ `roi` | `patient` and `feature` ∈ `mean_expression` | `celltype_frequency` (Section 2.1).

---

## 4. Storage Layer and Registry

Storage backends, registry schema, BioData hierarchy, Subject/Sample model, and MCP servers
are documented in full in `docs/registry.md` (symlinked as [[registry]] in the wiki).

Quick reference:

| Topic | Where |
|-------|-------|
| URI schemes, storage functions, data placement tiers | [[registry]] §1 |
| Registry schema tables, SQLite/PG setup | [[registry]] §2 |
| BioData JTI hierarchy and platform taxonomy | [[registry]] §2 + [[biodata-hierarchy]] |
| Subject/Sample model, 775-patient cohort schema | [[registry]] §2 + [[clinical-data-schema]] |
| MCP server tools (`sc-tools`, `sc-registry`) | [[registry]] §3 |

---

## 5. sc_tools vs Project-Specific (summary)

| Belongs to | Examples |
|------------|----------|
| **sc_tools** (generic) | `sc_tools.pl`, `sc_tools.tl`, `sc_tools.qc`, `sc_tools.pp`, `sc_tools.ingest`, `sc_tools.validate`, `sc_tools.bm`, `sc_tools.memory`, `sc_tools.utils` |
| **Project-specific** | `metadata/`, `results/`, `figures/`, `data/`, `outputs/`, project scripts |

---

## 6. Script Sanity Check

Scripts should use `$(PROJECT)` or equivalent for paths. Example: `$(PROJECT)/metadata/gene_signatures.json`, not `metadata/gene_signatures.json`.

### Active (in Snakemake dependency chain)

Scripts should write standard checkpoint names (Section 2.1). Legacy projects may still use old names during transition.

| Script | Phase slug | Primary output (standard) |
|--------|------------|----------------------------|
| Snakemake `spaceranger_count` rule | `ingest_raw` | `$(PROJECT)/data/{sample_id}/outs/` |
| `scripts/ingest.py` (per-sample loader) | `ingest_load` | `$(PROJECT)/data/{sample_id}/adata.ingested.h5ad` |
| `scripts/ingest.py` (concat + QC) | `qc_filter` | `$(PROJECT)/results/adata.filtered.h5ad` |
| `scripts/run_qc_report.py` | `qc_filter` | `$(PROJECT)/figures/QC/pre_filter_qc_YYYYMMDD.html` |
| Metadata join script | `metadata_attach` | `$(PROJECT)/results/adata.annotated.h5ad` |
| Preprocessing, clustering | `preprocess` | `$(PROJECT)/results/adata.normalized.h5ad` |
| Project scoring scripts (e.g. `projects/*/scripts/score_gene_signatures.py`) | `scoring` | `$(PROJECT)/results/adata.scored.h5ad` |
| `scripts/plot_deconvolution_spatial.py` | `scoring` | `$(PROJECT)/figures/deconvolution/` |
| `scripts/signature_heatmap_versioned.py` | `biology` | `$(PROJECT)/figures/manuscript/` |
| Manual cell typing workflow | `celltype_manual` | `$(PROJECT)/results/adata.celltyped.h5ad` |
| Downstream biology scripts | `biology` | `$(PROJECT)/figures/manuscript/` |
| Aggregation scripts | `meta_analysis` | `$(PROJECT)/results/adata.{level}.{feature}.h5ad` |

### Maintenance / utility scripts (active but not in a pipeline phase)

| Script | Purpose | Output |
|--------|---------|--------|
| `scripts/sync_wiki.py` | Regenerates `*.gen.md` wiki files from registry + pipeline DAG | `docs/wiki/**/*.gen.md` |
| `scripts/save_plan.py` | Hook called by ExitPlanMode to persist plans | `docs/wiki/plans/YYYY-MM-DD-{slug}.md` |
| `scripts/validate_checkpoint.py` | CLI wrapper around `sc_tools.validate`; exit 0/1 | Sentinel file or exit code |
| `scripts/run_qc_reports_all_projects.py` | Batch QC report generation across all registered projects | Per-project HTML reports |

### Legacy (read-only)
- **`scripts/old_code/`**: Reference only.

---

## 7. Operational Rules

1. **File placement:** All project outputs under `projects/<platform>/<project_name>/` (figures, results, metadata, data, outputs). No root-level metadata, results, or figures.
2. **Checkpoint nomenclature:** New pipelines and scripts **must** write AnnData checkpoints using the standard names in Section 2.1 (`adata.ingested.h5ad`, `adata.filtered.h5ad`, `adata.annotated.h5ad`, `adata.normalized.h5ad`, `adata.scored.h5ad`, `adata.celltyped.h5ad`). Legacy names (`adata.p0.h5ad`, `adata.raw.h5ad`, `adata.raw.p1.h5ad`, etc.) are accepted during transition but deprecated. Validators should check required metadata (Section 2.2) when loading checkpoints.
3. **Paths:** Scripts use `$(PROJECT)` or `PROJECT` variable for project paths.
4. **Statistics:** Benjamini–Hochberg (FDR); significance bars per `docs/skills.md`.
5. **Legacy:** Do not modify `scripts/old_code/`. Legacy checkpoint names (e.g. `adata.annotation.masked.h5ad`, `scvi.leiden.phenotyped.h5ad`) are allowed during migration; document which convention each project uses.
6. **Documentation:** Avoid apostrophes in generated text.
7. **Entry points:** Preprocessed projects may start at `preprocess`, `scoring`, or `celltype_manual`.

---

## 8. Testing

| Location | Scope | Fixtures |
|----------|-------|----------|
| `sc_tools/tests/` | Unit tests for sc_tools (pl, tl, qc, etc.) | Synthetic AnnData, in-memory |
| `projects/<platform>/<project>/tests/` | Integration/smoke tests for pipeline | Project fixtures, minimal data |

**Run tests:**
```bash
pytest sc_tools/tests/ -v
pytest projects/visium/ggo_visium/tests/ -v
# Or both: pytest sc_tools/tests/ projects/visium/ggo_visium/tests/ -v
```

**Implementation order:** (1) ggo_visium project tests, (2) sc_tools unit tests, (3) implement functions. All new code must compile and pass tests.

---

## 9. Development Environment

- **Environment:** `environment.yml` (conda env sc_tools), `requirements.txt` (pip). Package: `pip install -e ".[deconvolution]"`. Docker: see [project_setup.md](project_setup.md).
- **Libraries:** scanpy, squidpy, anndata, scvi-tools, tangram-sc; statannotations, pinguoin for figures.
