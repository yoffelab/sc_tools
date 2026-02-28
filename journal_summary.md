# Journal Summary: sc_tools Repository

Condensed summary of root `Journal.md` for context. Full entries and rationale are in `Journal.md`.

## Repo scope

- **sc_tools** is a reusable spatial omics toolkit and pipeline. Project-specific work lives under `projects/<data_type>/<project_name>/`. Root Journal covers repo structure, phasing, checkpoint names, and script roles only.

## Recent phase (2026-02)

- **Gene set scoring redesign (2026-02-28):** Full overhaul of Phase 3.5b gene set infrastructure. (1) `sc_tools/data/hallmark_human.json`: 50 bundled MSigDB Hallmark sets, stripped from robin metadata, offline. (2) `sc_tools/tl/gene_sets.py`: `load_hallmark`, `load_msigdb_json`, `load_gmt`, `list_gene_sets`, `validate_gene_signatures`, `merge_gene_signatures`, `update_gene_symbols`, `save_gene_signatures`. (3) `score_signature(method=...)`: `"scanpy"` (default), `"ucell"` (pyucell soft dep, rank AUC), `"ssgsea"` (gseapy soft dep); records `uns["scoring_method"]`. (4) `sc_tools/tl/gsea.py`: `run_ora` (Fisher exact + BH), `run_gsea_pseudobulk` (prerank via gseapy). (5) `sc_tools/pl/gsea.py`: `plot_gsea_dotplot`. Optional deps: `pip install sc-tools[geneset]`. All 80 tests pass; lint clean.
- **Generic deconvolution module (2026-02-26/27):** `sc_tools.tl.deconvolution()` with backend registry (cell2location, tangram, destvi). Per-library backed loading with memory profiling. Output: `obsm['cell_type_proportions']`, `obs['{method}_argmax']`. Fixed NaN bug for Seurat SCTransform references (negative values triggered double normalisation). Method-specific output naming: `adata.deconvolution.{method}.h5ad`. ggo_visium: both Tangram and Cell2location completed (29952 spots x 31 types). robin: Tangram done (280K spots x 39 types), Cell2location running on CPU.
- **Deconvolution spatial plots (2026-02-27):** `scripts/plot_deconvolution_spatial.py` saves per-library PNGs at 300 DPI alongside PDF. Method-specific output dirs (`figures/deconvolution/{method}/`). Spot sizing `120000/n_spots`; 98th percentile vmax for contrast; white background. Both Snakefiles updated (Phase 3.5b rules).
- **sc_tools skills as Cursor/Claude Code skill:** Repository `skills.md` exposed as `.cursor/skills/sc-tools-skills/SKILL.md` and `.claude/skills/sc-tools-skills/SKILL.md`; agent follows it for analysis and coding. Sandbox/local defaults: **Apptainer + Snakemake** (project_setup.md; skills.md §11). Docker used only for building the SIF image.
- **Container setup (current):** Docker builds image → `apptainer build containers/sc_tools.sif docker-daemon://sc_tools:latest`. All pipeline runs use `run_apptainer.sh` or Snakemake with `--use-apptainer`. On macOS, Apptainer requires a Linux VM or HPC; Docker is the local fallback.
- **Claude Code config added (2026-02-26):** CLAUDE.md at repo root; `.claude/skills/journal-and-mission/` and `.claude/skills/sc-tools-skills/` ported from Cursor. `.claude/settings.json` allowlists make, ruff, pytest, snakemake, git, apptainer.
- **Workflow:** journal_summary.md added everywhere; Mission.md is the todo list; in work mode the agent updates Mission after each prompt. Skill and Cursor rule enforce this; create_project.sh creates journal_summary.md for new projects.
- **CI/CD:** Linting (1) done. Snakemake (2) already adopted (each project has Snakefile). Order: 3 Sphinx, 4 PyPI, 5 GitHub Actions (last; testing complex). Containers: Apptainer (Linux/HPC, primary); Docker (macOS/Windows, fallback); auto-configure via `run_container.sh`.

- **Checkpoint nomenclature:** Standard filenames and required metadata are in Architecture.md Section 2 (e.g. `adata.raw.p1.h5ad`, `adata.normalized.p3.h5ad`, `adata.normalized.scored.p35.h5ad`). ggo_visium checkpoint files were renamed and scripts/Makefile updated to match.
- **Phasing:** Pipeline is 7 phases (1–7) with Phase 3.5 (Demographics) and 3.5b (Gene scoring, automated cell typing, deconvolution) branching from Phase 3. Phase 4 (Manual Cell Typing) is skippable if 3.5b is adequate. Gene signature storage: p35 uses `obsm['signature_score']` and `obsm['signature_score_z']` (full-path columns); `get_signature_df`/`get_signature_columns` read from obsm first, fall back to obs.
- **Layout:** All metadata, results, figures are project-specific. No root-level metadata/results/figures. Legacy scripts in `scripts/old_code/` per project. Single toolkit package: `sc_tools` only (ggo_tools/ggo_analysis removed).

## Key conventions

- Mission/Journal split: root = toolkit and pipeline; each project has its own Mission.md and Journal.md.
- New project: `./projects/create_project.sh <project_name> <data_type>`.
- Testing: two layers — package tests (`sc_tools/tests/`), project tests (`projects/.../tests/`). Order: ggo_visium tests first, then sc_tools tests, then implement.
