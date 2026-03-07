# sc_tools — Claude Code Configuration

## Sync Before Major Work

Read in order before any significant task:
1. @Mission.md — todo list and current priority
2. @journal_summary.md — condensed long-term memory
3. @Architecture.md — data flow, checkpoint names, directory layout

For subprojects, also read the project-level `Mission.md` and `journal_summary.md` under `projects/<platform>/<project>/`.

---

## Interaction Rules

### Work Mode (implementing, running code, making changes)
After EVERY prompt in work mode:
- Re-read Mission.md; mark completed steps `[x]`; update "In Progress"; add blockers/next steps

### Journal Workflow
- **Journal.md** — append dated entries after significant work (action, rationale, decisions)
- **journal_summary.md** — update condensed summary when Journal gets a new entry; keep it short
- Scope: repo root and each project with its own Journal.md

Plan mode (discussing only): updating Mission optional until execution begins.

---

## Container / Runtime

**Priority: Apptainer first (auto via `scripts/run_container.sh`):**
- Linux/HPC → Apptainer (primary)
- macOS/Windows → Docker (Apptainer not native; fallback)

**Running scripts:**
```bash
./scripts/run_container.sh projects/visium/ggo_visium python scripts/foo.py
./scripts/run_container.sh projects/visium/ggo_visium   # interactive bash
```

**Override:** `SC_TOOLS_RUNTIME=docker|apptainer ./scripts/run_container.sh ...`

**Build SIF:** `docker build -t sc_tools:latest . && apptainer build containers/sc_tools.sif docker-daemon://sc_tools:latest`

**Snakemake:** `snakemake -d projects/visium/ggo_visium -s projects/visium/ggo_visium/Snakefile <target>`

`run_apptainer.sh` and `run_docker.sh` are now thin wrappers around `run_container.sh`.

---

## Pipeline Phases (non-linear)

Phases are defined as semantic slugs in `sc_tools/pipeline.py`. Use `get_available_next(completed)` and `get_phase_checkpoint(slug)` to query the DAG programmatically.

| Slug | Old code | Name | Checkpoint | Required Data | QC Report |
|------|----------|------|------------|---------------|-----------|
| `ingest_raw` | 0a | Platform tools (Space Ranger / Xenium / IMC) | `data/{sample_id}/outs/` | Platform-specific raw output | |
| `ingest_load` | 0b | Load per-sample into AnnData / SpatialData | `data/{sample_id}/adata.h5ad` | `obs[sample, library_id, raw_data_dir]`, `obsm[spatial]`, `X` raw counts | |
| `qc_filter` | 1 | QC and Concatenation | `results/adata.raw.h5ad` | `obs[sample, raw_data_dir]`, `obsm[spatial]`, `X` raw, concatenated | `pre_filter_qc.html` |
| `metadata_attach` | 2 | Metadata Attachment (HIL) | `results/adata.annotated.h5ad` | + clinical columns in `obs` | `post_filter_qc.html` |
| `preprocess` | 3 | Preprocessing (+ integration benchmark) | `results/adata.normalized.h5ad` | `obsm[X_scvi or embedding]`, `obs[leiden]`, `adata.raw` backed up | `post_integration_qc.html` |
| `demographics` | 3.5 | Demographics (parallel branch, optional) | Figure 1 | Cohort metadata from `preprocess` | |
| `scoring` | 3.5b | Gene Scoring / Auto Cell Typing / Deconvolution | `results/adata.scored.h5ad` | `obsm[signature_score, signature_score_z]`, `uns[signature_score_report]` | |
| `celltype_manual` | 4 | Manual Cell Typing (HIL, optional, iterative) | `results/adata.celltyped.h5ad` | + `obs[celltype, celltype_broad]` | `post_celltyping_qc.html` |
| `biology` | 5 | Downstream Biology | `figures/manuscript/` | Reads from `scoring` or `celltype_manual` | |
| `meta_analysis` | 6/7 | Meta Analysis (optional) | `results/adata.{level}.{feature}.h5ad` | `obs` indexed by roi/patient; `X` = aggregated feature | |

All QC reports are date-versioned (`YYYYMMDD`) and saved to `figures/QC/`. Full validation contracts: Architecture.md Section 2.2.

Legacy p1/p2/p3/p35/p4 filenames still accepted during transition.

Entry points: `preprocess` (preprocessed), `scoring` (clustered), or `celltype_manual` (phenotyped). Full details: @Architecture.md §2.

---

## Testing Order

1. `projects/visium/ggo_visium/tests/` — project integration tests first
2. `sc_tools/tests/` — package unit tests second
3. Implement functions after tests pass

```bash
pytest sc_tools/tests/ -v
pytest projects/visium/ggo_visium/tests/ -v
make lint   # required before every commit
```

---

## Key Conventions

- **Paths:** All outputs under `projects/<platform>/<project>/` — no root-level results/figures/metadata
- **Checkpoint names:** Standard names only (see Architecture.md §2); no ad-hoc filenames
- **Script paths:** Use `$PROJECT_DIR` (set by run_container.sh) or Makefile `$(PROJECT)`
- **Signature scores:** `obsm['signature_score']` / `obsm['signature_score_z']`; not in `obs` by default
- **Colors:** `adata.uns[f'{name}_colors']`; never overwrite if length matches
- **No files > 1MB** in git
- **No apostrophes** in generated documentation text
- **Legacy:** `scripts/old_code/` read-only; refactor into sc_tools or new scripts
- **New project:** `./projects/create_project.sh <name> <data_type>`

---

## Analysis Standards

Full standards: @skills.md. Quick reference:
- FDR: Benjamini–Hochberg always; report adjusted p-values
- Significance bars: `*` < 0.05, `**` < 0.01, `***` < 0.001, `****` < 0.0001
- Group comparisons: 2 groups = pairwise; >2 groups = 1-vs-rest (default), all-pairwise as option
- Deconvolution: batch per `library_id` to avoid OOM
- Linting: Ruff (`make lint`) — never commit failing lint
- **Publication figures** (skills.md §12): 300+ DPI, Helvetica/Arial 5-8pt, color-blind safe palettes (Okabe-Ito default, no jet/rainbow), bold panel labels (lowercase for Nature, uppercase for Cell/Science), scale bars on spatial plots, exact P-values preferred. Use `marsilea` for complex composite figures.

---

## Skills

- `/journal-and-mission` — sync state, record decisions, update Mission.md
- `/sc-tools-skills` — apply analysis and coding standards from skills.md
