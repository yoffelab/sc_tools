# Projects

All analysis projects live under **data-type** folders: `visium/`, `visium_hd/`, `xenium/`, `imc/`.

Each project directory (e.g. `visium/ggo_visium/`) contains:

- `data/` — Raw sequencing and imaging
- `figures/` — Output visualizations (includes `QC/raw/`, `QC/post/`, `manuscript/`)
- `metadata/` — Gene signatures (JSON), sample metadata, `sample_metadata.csv`/`.xlsx`, `celltype_map.json`
- `scripts/` — Analysis scripts (and optionally `old_code/` for legacy)
- `results/` — Processed AnnData (.h5ad), CSVs
- `outputs/` — Intermediate tool outputs (e.g. deconvolution logs)
- `tests/` — Project integration tests (pytest)
- `Mission.md` — Project-specific goals (study aims, tasks, blockers)
- `Journal.md` — Project-specific analysis decision log (parameters, rationale, fixes)

Root `Mission.md` and `Journal.md` are for the toolkit and repo structure; project-specific aims (e.g. TLS, macrophage analysis for GGO Visium) live in each project's Mission and Journal.

---

## Pipeline Phases (Non-Linear)

The workflow has **human-in-the-loop** steps. See `README.md` (Pipeline Workflow section) for the full diagram.

| Phase | Name | Human-in-Loop? |
|-------|------|----------------|
| 1 | Data Ingestion & QC | No |
| 2 | Metadata Attachment | Yes (unless `sample_metadata.csv`/`.xlsx` provided) |
| 3 | Preprocessing | No |
| 3.5 | Demographics (branch) | Project-specific |
| 4 | Manual Cell Typing | Yes (iterative with `celltype_map.json`) |
| 5 | Downstream Biology | No |
| 6–7 | Meta Analysis (optional) | No |

**Entry points:** Preprocessed projects may start at Phase 3; already clustered projects may start at Phase 4.

---

## Testing

Each project can have `tests/` for integration/smoke tests. Run with `pytest projects/visium/ggo_visium/tests/ -v`. Implementation order: (1) project tests, (2) sc_tools package tests, (3) implement functions. See root `Mission.md` Section 5.

---

## Create a New Project

From repo root:

```bash
./create_project.sh <project_name> <data_type>
```

Example: `./create_project.sh ggo_visium visium` → `projects/visium/ggo_visium/` with the six subdirs.

**Migrate existing root-level content into ggo_visium (one-time):**

```bash
./migrate_to_ggo_visium.sh
```

**Run the pipeline:** From repo root, `make` uses `PROJECT=projects/visium/ggo_visium` by default. Override with `make PROJECT=projects/visium/other_project phase1`, etc.
