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

| Phase | Name | Checkpoint |
|-------|------|------------|
| 1 | Data Ingestion & QC | `adata.raw.p1.h5ad` |
| 2 | Metadata Attachment (HIL) | `adata.annotated.p2.h5ad` |
| 3 | Preprocessing | `adata.normalized.p3.h5ad` |
| 3.5 | Demographics (parallel branch) | Figure 1 |
| 3.5b | Gene Scoring / Auto Cell Typing / Deconvolution | `adata.normalized.scored.p35.h5ad` |
| 4 | Manual Cell Typing (HIL, skippable) | `adata.celltyped.p4.h5ad` |
| 5 | Downstream Biology | `figures/manuscript/` |
| 6–7 | Meta Analysis (optional) | `adata.{level}.{feature}.h5ad` |

Entry points: Phase 3 (preprocessed), Phase 3.5b (clustered), or Phase 4 (phenotyped). Full details: @Architecture.md §2.

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

---

## Skills

- `/journal-and-mission` — sync state, record decisions, update Mission.md
- `/sc-tools-skills` — apply analysis and coding standards from skills.md
