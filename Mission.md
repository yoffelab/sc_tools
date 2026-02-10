# Mission: sc_tools Toolkit and Pipeline

**Scope:** Repository-level and generalizable goals. Project-specific objectives (e.g. TLS, macrophage analysis) live in each project's `Mission.md` under `projects/<data_type>/<project_name>/Mission.md`.

**Last Updated:** 2025-02-09

---

## 1. Toolkit and Pipeline (General)

- **Pipeline phases (1–7):** Non-linear workflow with human-in-loop steps. See `README.md` for the diagram and phase summary.
- **sc_tools:** Reusable package (`sc_tools.pl`, `sc_tools.tl`, `sc_tools.qc`, etc.) for QC, plotting, testing, colocalization, I/O. Keep generic; project-specific logic stays in project scripts.
- **Reproducibility:** Makefile is project-aware (`PROJECT ?= projects/visium/ggo_visium`). Each project has `data/`, `figures/`, `metadata/`, `scripts/`, `results/`, `outputs/` under `projects/<platform>/<project_name>/`.
- **Standards:** All projects follow `skills.md` (statistics, significance bars, FDR). Documentation avoids apostrophes.

---

## 2. General Analysis — To-Do and Implementation Status

All paths below are project-specific: `projects/<platform>/<project_name>/...`. Nothing lives at repo root for `metadata/`, `results/`, or `figures/`.

### Phase 1: Data Ingestion & QC

#### Ingestion (Platform-Specific)
- [ ] **Visium / Visium HD:** fastq + H&E images + Cytassist images → Space Ranger → cloupe. Convert cloupe to AnnData. Keep H&E images when concatenating.
- [ ] **IMC:** mcd, txt files → segmentation and concatenation → h5ad.
- [ ] **Xenium:** Assume preprocessed; load into AnnData.
- [ ] **Required annotations:** `adata.obs['sample']` (sample of origin), `adata.obs['raw_data_dir']` (backup path for original data). Spatial coords in `adata.obsm['spatial']`.

#### QC Metrics (sc_tools.qc)
- [ ] Implement `calculate_qc_metrics` (wrap scanpy `pp.calculate_qc_metrics`).
- [ ] Implement `filter_cells`, `filter_genes` (wrap scanpy).
- [ ] Implement `highly_variable_genes` (wrap scanpy).
- [ ] Implement `spatially_variable_genes` (wrap squidpy). Sample-specific QC on spots/cells and genes/proteins.
- [ ] Two QC versions: pre-normalization and post-normalization.

#### QC Plots
- [ ] 2x2 grid: total_count (total expression), gene/protein expression histogram, log1p version (all platforms).
- [ ] Spatial transcriptomics: % mt, % hb per spot.
- [ ] Multipage spatial plot per sample: total_count, log1p, %mt (1x3 subplot) as QC report.
- [ ] Save under `projects/<platform>/<project>/figures/QC/raw/` (pre) and `figures/QC/post/` (post).
- [ ] **Future:** MA plots (per gene/protein across samples).

---

### Phase 2: Metadata Attachment (Human-in-Loop)

- [ ] **Bypass file:** `projects/<platform>/<project>/metadata/sample_metadata.csv` or `.xlsx` mapping `sample` → clinical columns. sc_tools guides the join when provided.
- [ ] Without file: human prepares map; cannot skip via automatic pipelining.
- [ ] Must be done as early as possible in the project cycle.
- [ ] Preprocessed projects may skip to Phase 3 or 4 if already annotated.

---

### Phase 3: Preprocessing

- [ ] **Before any modification:** Backup `adata.raw`.
- [ ] Filter cells, genes, proteins, samples failing QC (project-specific thresholds).
- [ ] Normalize, batch correct (e.g. scVI), run automated clustering.
- [ ] Post-normalization QC report → `figures/QC/post/`.
- [ ] **No automated cell typing here** (moved to Phase 3.5b).

---

### Phase 3.5: Demographics (Branching)

- [ ] **sc_tools helpers:** piechart, histogram, violinplot, kdeplot, barplot, stacked barplot, scatterplot, correlogram, heatmap.
- [ ] Figure 1 for population-based studies describing the cohort.

---

### Phase 3.5b: Gene Scoring, Automated Cell Typing, Deconvolution

- [ ] **Gene signature storage (TODO):** Store scores in **`adata.obsm['sig:hallmark']`** for Hallmark processes; for each project-specific `metadata/{signature_name}.json`, store in **`adata.obsm['sig:{signature_name}']`**. Do not default to storing in `adata.obs`; use `obsm` for signature matrices.
- [ ] **Always apply:** Basic gene sets (e.g. Hallmark) plus any signatures provided in project `metadata/*.json`.
- [ ] Automated cell typing (cluster → celltype); non-transformer and transformer models.
- [ ] Optional cell-type deconvolution (Tangram, Cell2location, DestVI; batch per library/sample).
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

- **GGO Visium:** `projects/visium/ggo_visium/Mission.md`.
- **Other projects:** Add `Mission.md` and `Journal.md` via `./projects/create_project.sh <project_name> <data_type>`.
- **Entry points:** Preprocessed projects may start at Phase 3 or 4.

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

- **sc_tools package:** `pl/` (spatial, heatmaps, statistical, volcano, save), `tl/` (testing, colocalization, deconvolution, io), `memory/` (profiling, gpu), `qc/` (placeholder).
- **projects layout:** `visium/`, `visium_hd/`, `xenium/`, `imc/`, `cosmx/`; each project has `data/`, `figures/`, `metadata/`, `scripts/`, `results/`, `outputs/`.
- **create_project.sh:** `./projects/create_project.sh <project_name> <data_type>`.
- **Makefile project-aware:** `PROJECT ?= projects/visium/ggo_visium`; paths use `$(PROJECT)/...`.
- **Metadata moved:** Root `metadata/` moved to `projects/visium/ggo_visium/metadata/`.
- **Legacy merged:** `scripts/old_code/` (read-only reference).

### In Progress — Next immediate steps

**Testing (order: 1st ggo_visium, 2nd sc_tools, 3rd functions)**
- [ ] ggo_visium project tests: `projects/visium/ggo_visium/tests/`. Integration/smoke tests for Makefile and scripts.
- [ ] sc_tools package tests: `sc_tools/tests/`. Unit tests for pl, tl, qc with synthetic fixtures.
- [ ] Add pytest to pyproject.toml / requirements.

**Makefile & scripts**
- [ ] Align Makefile targets with Phases 1–7 (see README Pipeline Workflow).
- [ ] Update scripts to use `$(PROJECT)/metadata/` instead of `metadata/`. Affected: score_gene_signatures.py, tumor_differences.py, tls_analysis.py, manuscript_spatial_plots.py, celltyping.py, etc.
- [ ] Modular scripts: config-driven; import from `sc_tools`; thin orchestration only.

### To Do (later)

- [ ] Organize production scripts by phase within project.
- [ ] Update imports to use `sc_tools.pl.*`, `sc_tools.tl.*` everywhere.
- [ ] imc-analysis: review, identify functionalities, integrate into `sc_tools`.
- [ ] Documentation: migration guide, docstrings, API docs.

### Operational notes

- **API:** `st.pl.*`, `st.tl.*`, `st.qc.*`, `st.memory.*` (scanpy-style).
- **New project:** `./projects/create_project.sh <project_name> visium|visium_hd|xenium|imc|cosmx`. Run make with `PROJECT=projects/<type>/<name>`.
- **Legacy:** `scripts/old_code/` is read-only. Do not modify; refactor into new scripts or `sc_tools` if needed.

---

## 7. Reference

- **README.md** — Pipeline workflow diagram, phase summary.
- **Architecture.md** — Directory layout, phase details, project-specific paths, testing structure.
- **skills.md** — Coding and statistical standards.
- **Per-project:** `projects/<data_type>/<project_name>/Mission.md` and `Journal.md`.
