# Project Architecture: Scalable Spatial Omics Analysis

This document outlines the directory structure, data flow, and script inventory. The layout is **scalable across projects and data types** (Visium, Visium HD, Xenium, IMC, CosMx). All **metadata**, **results**, and **figures** are **project-specific** under `projects/<platform>/<project_name>/`. There is no repo-root `metadata/`, `results/`, or `figures/`. The pipeline is non-linear; see `README.md` (Pipeline Workflow section) for the diagram.

## 1. Directory Overview

```text
.
├── Architecture.md         # System roadmap (this file)
├── Mission.md              # Todo list and roadmap (general); project-specific in projects/<type>/<name>/Mission.md
├── Journal.md              # Repo-level decision log; project-specific in projects/<type>/<name>/Journal.md
├── journal_summary.md      # Short summary of Journal.md for context; per-project under projects/<type>/<name>/
├── skills.md               # Mandatory coding and statistical standards
├── pyproject.toml          # Package build (sc_tools installable)
├── environment.yml         # Conda environment
├── requirements.txt        # pip dependencies
├── Makefile                # Pipeline orchestration (project-aware; default PROJECT=projects/visium/ggo_visium)
├── ENVIRONMENT_SETUP.md    # Environment setup notes
│
├── sc_tools/               # Reusable Python package (scanpy-style API) — NOT project-specific
│   ├── pl/                  # Plotting: spatial, heatmaps, statistical, volcano, save
│   ├── tl/                  # Tools: testing, colocalization, deconvolution, io
│   ├── qc/                  # QC: calculate_qc_metrics, filter_cells, filter_genes, highly_variable_genes, spatially_variable_genes
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
│   ├── xenium/
│   ├── imc/
│   └── cosmx/
└── .gitignore, .githooks/
```

**Creating a new project:** Run `./projects/create_project.sh <project_name> <data_type>`. Valid `data_type`: `visium` | `visium_hd` | `xenium` | `imc` | `cosmx`.

---

## 2. Project-Specific Paths and Checkpoint Nomenclature

All of the following live under `projects/<platform>/<project_name>/`. **Checkpoint AnnData files must use the standard names below** so tooling can validate metadata and downstream scripts can assume consistent paths.

### 2.1 Mandatory checkpoint filenames (results/)

| Phase | Standard path | Description |
|-------|----------------|-------------|
| **1** | `results/adata.raw.p1.h5ad` | Raw unnormalized AnnData after ingestion and QC. |
| **2** | `results/adata.annotated.p2.h5ad` | AnnData with clinical/sample metadata (and masks if applicable) attached. |
| **3** | `results/adata.normalized.p3.h5ad` | Filtered, normalized, batch-corrected, clustered (no cell typing). |
| **3.5b** | `results/adata.normalized.scored.p35.h5ad` | Normalized + gene signature scores (and automated cell typing if applied). |
| **4** | `results/adata.celltyped.p4.h5ad` | Manual cell type labels applied (celltype, celltype_broad). |
| **5–7** (aggregated) | `results/adata.{level}.{feature}.h5ad` | Aggregated by `level` ∈ `roi` | `patient` and `feature` ∈ `mean_expression` | `celltype_frequency`. |

**Aggregated examples:** `results/adata.roi.mean_expression.h5ad`, `results/adata.patient.celltype_frequency.h5ad`.

### 2.2 Required metadata (for validation)

Checkpoints must satisfy the following so scripts and validators can rely on them. Tooling (e.g. `sc_tools` or a project test) should **check** these before using a checkpoint.

| Checkpoint | Required `adata` contents |
|------------|---------------------------|
| **adata.raw.p1.h5ad** | `obs['sample']`, `obs['raw_data_dir']` (or equivalent), `obsm['spatial']`; `X` raw counts; no normalization. |
| **adata.annotated.p2.h5ad** | All of p1; `obs` includes clinical/metadata columns from `metadata/sample_metadata.csv` (or join equivalent); optional `uns` keys for images/masks. |
| **adata.normalized.p3.h5ad** | Normalized/batch-corrected representation (e.g. `obsm['X_scvi']`); `obs['leiden']` (or cluster column); `adata.raw` backed up. |
| **adata.normalized.scored.p35.h5ad** | All of p3; signature scores in `obsm['sig:hallmark']` and/or `obsm['sig:{name}']` per project `metadata/{name}.json` (target); optional `obs['celltype']`/`celltype_broad` if automated typing run. |
| **adata.celltyped.p4.h5ad** | All of p35; `obs['celltype']`, `obs['celltype_broad']` from `metadata/celltype_map.json`. |
| **adata.{level}.{feature}.h5ad** | `obs` indexed by `level` (roi or patient); `X` or layer holds aggregated `feature` (mean expression or celltype frequency). |

### 2.3 Other project paths (metadata, optional results, figures)

| Path | Description |
|------|-------------|
| `metadata/sample_metadata.csv` or `.xlsx` | Sample→clinical metadata map. Enables Phase 2 bypass. |
| `metadata/celltype_map.json` | cluster_id→celltype mapping for Phase 4. |
| `metadata/gene_signatures.json` | Gene signatures for scoring (and per-signature `metadata/{name}.json` for obsm storage). |
| `results/adata.deconvolution.h5ad` | Cell-type proportions (optional; Phase 3.5b). |
| `figures/QC/raw/` | Pre-normalization QC reports. |
| `figures/QC/post/` | Post-normalization QC reports. |
| `figures/manuscript/` | Publication figures. |

**Legacy / migration:** Existing projects may still use `adata.annotation.masked.h5ad`, `scvi.leiden.phenotyped.h5ad`, `adata.img.genescores.h5ad` until scripts are updated. New pipelines and scripts **must** write the standard checkpoint names above. Validators should accept either legacy or standard names and report which convention is used.

**Makefile:** Use `$(PROJECT)` for all project paths, e.g. `$(PROJECT)/results/adata.raw.p1.h5ad`.

---

## 3. Phase Details and Data Flow

### Phase 1: Data Ingestion & QC

**Platform-specific ingestion:**
- **Visium / Visium HD:** fastq, H&E, Cytassist → Space Ranger → cloupe → AnnData. Keep H&E when concatenating.
- **IMC:** mcd, txt → segmentation + concatenation → h5ad.
- **Xenium:** Assume preprocessed.

**Required annotations:** `adata.obs['sample']`, `adata.obs['raw_data_dir']`, `adata.obsm['spatial']`.

**QC (sc_tools.qc):** scanpy `calculate_qc_metrics`, `filter_cells`, `filter_genes`, `highly_variable_genes`; squidpy spatially variable genes. 2x2 grid, % mt/% hb for spatial, multipage spatial (total_count, log1p, %mt). Pre → `figures/QC/raw/`; post → `figures/QC/post/`.

**Outputs:** `$(PROJECT)/results/adata.raw.p1.h5ad`, `$(PROJECT)/figures/QC/raw/*`. Must satisfy required metadata for p1 (Section 2.2).

---

### Phase 2: Metadata Attachment (Human-in-Loop)

**Bypass:** Provide `$(PROJECT)/metadata/sample_metadata.csv` or `.xlsx`.

**Without file:** Human prepares map; cannot skip automatically.

**Outputs:** `$(PROJECT)/results/adata.annotated.p2.h5ad`. Must satisfy required metadata for p2 (Section 2.2).

---

### Phase 3: Preprocessing

Backup `adata.raw`; filter; normalize; batch correct; cluster. Post-QC → `$(PROJECT)/figures/QC/post/`. No automated cell typing (that is in Phase 3.5b).

**Outputs:** `$(PROJECT)/results/adata.normalized.p3.h5ad`. Must satisfy required metadata for p3 (Section 2.2).

---

### Phase 3.5: Demographics (Branching)

Separate branch from Phase 3 (parallel to 3.5b). sc_tools helpers: piechart, histogram, violinplot, barplot, stacked barplot, scatterplot, correlogram, heatmap. Figure 1 for cohort description.

---

### Phase 3.5b: Gene Scoring, Automated Cell Typing, Deconvolution

Separate branch from Phase 3 (parallel to 3.5); connects to Phase 4. Always apply basic gene sets (e.g. Hallmark) and any project-provided signatures from `metadata/{signature_name}.json`. Store scores in **`adata.obsm['sig:hallmark']`**, **`adata.obsm['sig:{signature_name}']`** (not in `obs` by default). Automated cell typing (cluster → celltype). Optional cell-type deconvolution (DestVI, Cell2location, Tangram) → `$(PROJECT)/results/adata.deconvolution.h5ad`.

**Outputs:** `$(PROJECT)/results/adata.normalized.scored.p35.h5ad` (must satisfy Section 2.2); optional `adata.deconvolution.h5ad`. Required for Phase 5.

---

### Phase 4: Manual Cell Typing (Human-in-Loop; Skippable)

Skippable if automated cell typing in 3.5b is adequate. JSON format `{cluster_id: {celltype_name, total_obs_count}}`; match cluster_id type; produce celltype and celltype_broad. Iterative until satisfactory. Save to `$(PROJECT)/metadata/celltype_map.json`.

**Outputs:** `$(PROJECT)/results/adata.celltyped.p4.h5ad`. Must satisfy required metadata for p4 (Section 2.2).

---

### Phase 5: Downstream Biology

Uses gene scores and (optionally) deconvolution from Phase 3.5b. Spatial/process analysis, colocalization, neighborhood enrichment, publication figures. Reads from `adata.normalized.scored.p35.h5ad` (or p4).

---

### Phase 6–7: Meta Analysis (Optional)

Aggregate ROI/patient; downstream on aggregated data.

**Outputs:** `$(PROJECT)/results/adata.{level}.{feature}.h5ad` with `level` ∈ `roi` | `patient` and `feature` ∈ `mean_expression` | `celltype_frequency` (Section 2.1).

---

## 4. sc_tools vs Project-Specific

| Belongs to | Examples |
|------------|----------|
| **sc_tools** (generic) | `sc_tools.pl`, `sc_tools.tl`, `sc_tools.qc`, `sc_tools.memory`, `sc_tools.utils` |
| **Project-specific** | `metadata/`, `results/`, `figures/`, `data/`, `outputs/`, project scripts |

---

## 5. Script Sanity Check

Scripts should use `$(PROJECT)` or equivalent for paths. Example: `$(PROJECT)/metadata/gene_signatures.json`, not `metadata/gene_signatures.json`.

### Active (in Makefile dependency chain)

Scripts should write standard checkpoint names (Section 2.1). Legacy projects may still use old names during migration.

| Script | Phase | Primary output (standard) |
|--------|-------|----------------------------|
| Platform-specific ingestion | 1 | $(PROJECT)/results/adata.raw.p1.h5ad |
| QC script (or sc_tools.qc) | 1 | $(PROJECT)/figures/QC/raw/ |
| Metadata join script | 2 | $(PROJECT)/results/adata.annotated.p2.h5ad |
| preprocessing, clustering | 3 | $(PROJECT)/results/adata.normalized.p3.h5ad |
| score_gene_signatures, automated celltyping, deconvolution | 3.5b | $(PROJECT)/results/adata.normalized.scored.p35.h5ad |
| Manual cell typing workflow | 4 | $(PROJECT)/results/adata.celltyped.p4.h5ad |
| tumor_differences, process_colocalization, etc. | 5 | $(PROJECT)/figures/manuscript/ |
| Aggregation scripts | 6–7 | $(PROJECT)/results/adata.{level}.{feature}.h5ad |

### Legacy (read-only)
- **`scripts/old_code/`**: Reference only.

---

## 6. Operational Rules

1. **File placement:** All project outputs under `projects/<platform>/<project_name>/` (figures, results, metadata, data, outputs). No root-level metadata, results, or figures.
2. **Checkpoint nomenclature:** New pipelines and scripts **must** write AnnData checkpoints using the standard names in Section 2.1 (`adata.raw.p1.h5ad`, `adata.annotated.p2.h5ad`, etc.). Validators should check required metadata (Section 2.2) when loading checkpoints.
3. **Paths:** Scripts use `$(PROJECT)` or `PROJECT` variable for project paths.
4. **Statistics:** Benjamini–Hochberg (FDR); significance bars per `skills.md`.
5. **Legacy:** Do not modify `scripts/old_code/`. Legacy checkpoint names (e.g. `adata.annotation.masked.h5ad`, `scvi.leiden.phenotyped.h5ad`) are allowed during migration; document which convention each project uses.
6. **Documentation:** Avoid apostrophes in generated text.
7. **Entry points:** Preprocessed projects may start at Phase 3, 3.5b, or 4.

---

## 7. Testing

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

## 8. Development Environment

- **Environment:** `environment.yml` (conda), `requirements.txt` (pip). Package: `pip install -e ".[deconvolution]"`.
- **Libraries:** scanpy, squidpy, anndata, scvi-tools, tangram-sc; statannotations, pinguoin for figures.
