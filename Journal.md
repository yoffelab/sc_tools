# Research Journal & Decision Log: sc_tools Repository

This journal documents **repository-level** technical and structural decisions. Project-specific analysis decisions (e.g. GGO Visium TLS, macrophage, tumor differences) are in each project's `Journal.md` under `projects/<data_type>/<project_name>/Journal.md`.

## Repo-wide notes
* **AI-assisted workflow:** Agentic scientific programming is used; project content lives under `projects/<data_type>/<project_name>/`.
* **Constraint:** Generated documentation and communications avoid apostrophes.

---

## Log Entries (toolkit / repo structure)

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
