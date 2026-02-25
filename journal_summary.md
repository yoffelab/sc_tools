# Journal Summary: sc_tools Repository

Condensed summary of root `Journal.md` for context. Full entries and rationale are in `Journal.md`.

## Repo scope

- **sc_tools** is a reusable spatial omics toolkit and pipeline. Project-specific work lives under `projects/<data_type>/<project_name>/`. Root Journal covers repo structure, phasing, checkpoint names, and script roles only.

## Recent phase (2026-02)

- **Docker + conda + UV:** Dockerfile uses continuumio/miniconda3 with conda env sc_tools; UV installs sc_tools via `--python` flag. project_setup.md documents build, run, per-project usage (ggo_visium, robin, lymph_dlbcl), local dev. environment.yml renamed to sc_tools; Robin has run_docker.sh.
- **Workflow:** journal_summary.md added everywhere; Mission.md is the todo list; in work mode the agent updates Mission after each prompt. Skill and Cursor rule enforce this; create_project.sh creates journal_summary.md for new projects.
- **CI/CD:** Linting (1) done. Order: 2 Nextflow (early for easier porting), 3 Sphinx, 4 PyPI, 5 GitHub Actions (last; testing complex). Nextflow + Docker (local) / Singularity (HPC); auto-configure based on environment.

- **Checkpoint nomenclature:** Standard filenames and required metadata are in Architecture.md Section 2 (e.g. `adata.raw.p1.h5ad`, `adata.normalized.p3.h5ad`, `adata.normalized.scored.p35.h5ad`). ggo_visium checkpoint files were renamed and scripts/Makefile updated to match.
- **Phasing:** Pipeline is 7 phases (1–7) with Phase 3.5 (Demographics) and 3.5b (Gene scoring, automated cell typing, deconvolution) branching from Phase 3. Phase 4 (Manual Cell Typing) is skippable if 3.5b is adequate. Gene signature storage: p35 uses `obsm['signature_score']` and `obsm['signature_score_z']` (full-path columns); `get_signature_df`/`get_signature_columns` read from obsm first, fall back to obs.
- **Layout:** All metadata, results, figures are project-specific. No root-level metadata/results/figures. Legacy scripts in `scripts/old_code/` per project. Single toolkit package: `sc_tools` only (ggo_tools/ggo_analysis removed).

## Key conventions

- Mission/Journal split: root = toolkit and pipeline; each project has its own Mission.md and Journal.md.
- New project: `./projects/create_project.sh <project_name> <data_type>`.
- Testing: two layers — package tests (`sc_tools/tests/`), project tests (`projects/.../tests/`). Order: ggo_visium tests first, then sc_tools tests, then implement.
