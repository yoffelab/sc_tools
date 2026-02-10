# Research Journal & Decision Log: sc_tools Repository

This journal documents **repository-level** technical and structural decisions. Project-specific analysis decisions (e.g. GGO Visium TLS, macrophage, tumor differences) are in each project's `Journal.md` under `projects/<data_type>/<project_name>/Journal.md`.

## Repo-wide notes
* **AI-assisted workflow:** Agentic scientific programming is used; project content lives under `projects/<data_type>/<project_name>/`.
* **Constraint:** Generated documentation and communications avoid apostrophes.

---

## Log Entries (toolkit / repo structure)

### [2025-02-09] - Pipeline Phasing Overhaul (Phases 1–7)
- **Action:** Rewrote pipeline into 7 phases with non-linear workflow and human-in-loop steps.
- **Phases:** (1) Data Ingestion & QC; (2) Metadata Attachment (HIL unless sample_metadata.csv); (3) Preprocessing; (3.5) Demographics; (4) Manual Cell Typing (HIL, iterative); (5) Downstream Biology; (6–7) Meta Analysis (optional).
- **Additions:** Created `WORKFLOW.md` with Mermaid diagram; `sc_tools.qc` placeholder (metrics, spatial, plots); project `metadata/sample_metadata.csv`, `metadata/celltype_map.json` as bypass files.

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
