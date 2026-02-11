# Research Journal & Decision Log: sc_tools Repository

This journal documents **repository-level** technical and structural decisions. Project-specific analysis decisions (e.g. GGO Visium TLS, macrophage, tumor differences) are in each project's `Journal.md` under `projects/<data_type>/<project_name>/Journal.md`.

## Repo-wide notes
* **AI-assisted workflow:** Agentic scientific programming is used; project content lives under `projects/<data_type>/<project_name>/`.
* **Constraint:** Generated documentation and communications avoid apostrophes.

---

## Log Entries (toolkit / repo structure)

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
