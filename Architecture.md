# Project Architecture: Scalable Spatial Omics Analysis

This document outlines the directory structure, data flow, and script inventory. The layout is **scalable across projects and data types** (Visium, Visium HD, Xenium, IMC). All **metadata**, **results**, and **figures** are **project-specific** under `projects/<platform>/<project_name>/`. There is no repo-root `metadata/`, `results/`, or `figures/`. The pipeline is non-linear; see `WORKFLOW.md` for the diagram.

## 1. Directory Overview

```text
.
├── Architecture.md         # System roadmap (this file)
├── Mission.md              # Toolkit and pipeline (general); project-specific in projects/<type>/<name>/Mission.md
├── WORKFLOW.md             # Non-linear pipeline diagram (Mermaid) and phase summary
├── Journal.md              # Repo-level decision log; project-specific in projects/<type>/<name>/Journal.md
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
│   │       └── Journal.md
│   ├── visium_hd/
│   ├── xenium/
│   └── imc/
└── .gitignore, .githooks/
```

**Creating a new project:** Run `./projects/create_project.sh <project_name> <data_type>`. Valid `data_type`: `visium` | `visium_hd` | `xenium` | `imc`.

---

## 2. Project-Specific Paths (Explicit)

All of the following live under `projects/<platform>/<project_name>/`:

| Path | Description |
|------|-------------|
| `metadata/sample_metadata.csv` or `.xlsx` | Sample→clinical metadata map. Enables Phase 2 bypass. |
| `metadata/celltype_map.json` | cluster_id→celltype mapping for Phase 4. |
| `metadata/gene_signatures.json` | Gene signatures for scoring. |
| `results/adata.raw.h5ad` | Raw unnormalized AnnData after Phase 1. |
| `results/adata.annotation.masked.h5ad` | AnnData with masks (project-specific). |
| `results/scvi.h5ad`, `results/scvi.leiden.h5ad` | Batch-corrected, clustered AnnData. |
| `results/scvi.leiden.phenotyped.h5ad` | Clustered + celltype annotated. |
| `results/adata.img.genescores.h5ad` | Gene signature scores. |
| `results/adata.deconvolution.h5ad` | Cell-type proportions (optional). |
| `figures/QC/raw/` | Pre-normalization QC reports. |
| `figures/QC/post/` | Post-normalization QC reports. |
| `figures/manuscript/` | Publication figures. |

**Makefile:** Use `$(PROJECT)` for all project paths, e.g. `$(PROJECT)/metadata/gene_signatures.json`, `$(PROJECT)/results/adata.raw.h5ad`.

---

## 3. Phase Details and Data Flow

### Phase 1: Data Ingestion & QC

**Platform-specific ingestion:**
- **Visium / Visium HD:** fastq, H&E, Cytassist → Space Ranger → cloupe → AnnData. Keep H&E when concatenating.
- **IMC:** mcd, txt → segmentation + concatenation → h5ad.
- **Xenium:** Assume preprocessed.

**Required annotations:** `adata.obs['sample']`, `adata.obs['raw_data_dir']`, `adata.obsm['spatial']`.

**QC (sc_tools.qc):** scanpy `calculate_qc_metrics`, `filter_cells`, `filter_genes`, `highly_variable_genes`; squidpy spatially variable genes. 2x2 grid, % mt/% hb for spatial, multipage spatial (total_count, log1p, %mt). Pre → `figures/QC/raw/`; post → `figures/QC/post/`.

**Outputs:** `$(PROJECT)/results/adata.raw.h5ad`, `$(PROJECT)/figures/QC/raw/*`.

---

### Phase 2: Metadata Attachment (Human-in-Loop)

**Bypass:** Provide `$(PROJECT)/metadata/sample_metadata.csv` or `.xlsx`.

**Without file:** Human prepares map; cannot skip automatically.

---

### Phase 3: Preprocessing

Backup `adata.raw`; filter; normalize; batch correct; cluster; automated cell typing. Post-QC → `$(PROJECT)/figures/QC/post/`.

---

### Phase 3.5: Demographics (Branching)

sc_tools helpers: piechart, histogram, violinplot, barplot, stacked barplot, scatterplot, correlogram, heatmap. Figure 1 for cohort description.

---

### Phase 4: Manual Cell Typing (Human-in-Loop)

JSON format `{cluster_id: {celltype_name, total_obs_count}}`; match cluster_id type; produce celltype and celltype_broad. Iterative until satisfactory. Save to `$(PROJECT)/metadata/celltype_map.json`.

---

### Phase 5: Downstream Biology

Gene scoring, deconvolution, spatial/process analysis, colocalization, neighborhood enrichment.

---

### Phase 6–7: Meta Analysis (Optional)

Aggregate ROI/patient; downstream on aggregated data.

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
| Script | Phase | Primary output |
|--------|-------|----------------|
| Platform-specific ingestion | 1 | $(PROJECT)/results/adata.raw.h5ad |
| QC script (or sc_tools.qc) | 1 | $(PROJECT)/figures/QC/raw/ |
| Metadata join script | 2 | AnnData with clinical metadata |
| preprocessing, clustering, celltyping | 3 | $(PROJECT)/results/scvi.leiden.phenotyped.h5ad |
| score_gene_signatures, deconvolution | 5 | $(PROJECT)/results/adata.img.genescores.h5ad |
| tumor_differences, process_colocalization, etc. | 5 | $(PROJECT)/figures/manuscript/ |
| Aggregation scripts | 6–7 | ROI/patient aggregated data |

### Legacy (read-only)
- **`scripts/old_code/`**: Reference only.

---

## 6. Operational Rules

1. **File placement:** All project outputs under `projects/<platform>/<project_name>/` (figures, results, metadata, data, outputs). No root-level metadata, results, or figures.
2. **Paths:** Scripts use `$(PROJECT)` or `PROJECT` variable for project paths.
3. **Statistics:** Benjamini–Hochberg (FDR); significance bars per `skills.md`.
4. **Legacy:** Do not modify `scripts/old_code/`.
5. **Documentation:** Avoid apostrophes in generated text.
6. **Entry points:** Preprocessed projects may start at Phase 3 or 4.

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
