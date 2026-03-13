# Projects

All analysis projects live under **data-type** folders: `visium/`, `visium_hd/`, `xenium/`, `imc/`, `cosmx_1k/`, `cosmx_6k/`, `cosmx_full_library/`.

## Create a New Project

From the repository root:

```bash
./projects/create_project.sh <project_name> <data_type>
```

Valid `data_type`: `visium` | `visium_hd` | `visium_hd_cell` | `xenium` | `imc` | `cosmx_1k` | `cosmx_6k` | `cosmx_full_library`

Example: `./projects/create_project.sh my_study visium_hd` creates `projects/visium_hd/my_study/` with the full directory tree.

---

## Project Directory Structure

Each project directory contains:

| Path | Purpose |
|------|---------|
| `data/` | Raw sequencing and imaging (or symlinks to HPC) |
| `metadata/` | Gene signatures (JSON), `sample_metadata.csv`/`.xlsx`, `celltype_map.json`, `phase0/` batch manifests |
| `results/` | AnnData checkpoints (.h5ad) |
| `figures/` | QC (`QC/raw/`, `QC/post/`) and manuscript figures |
| `scripts/` | Project-specific analysis scripts |
| `outputs/` | Intermediate tool outputs (deconvolution logs, etc.) |
| `tests/` | Project integration tests (pytest) |
| `Snakefile` | Project pipeline |
| `config.yaml` | Snakemake config (paths, parameters) |
| `Plan.md` | Master todo list, phase status, and implementation plan |
| `Journal.md` | Project-specific decision log |

`docs/Plan.md` and `docs/Journal.md` cover toolkit-level work; each project has its own for study-specific goals.

---

## Pipeline Phases

| Phase | Name | Human-in-Loop? |
|-------|------|----------------|
| **0a** | Platform tools (Space Ranger / Xenium Ranger / IMC) | No |
| **0b** | Load per-sample into AnnData | No |
| **1** | QC and Concatenation | No |
| **2** | Metadata Attachment | Yes (unless `sample_metadata.csv`/`.xlsx` provided) |
| **3** | Preprocessing | No |
| **3.5** | Demographics (parallel branch) | Project-specific |
| **3.5b** | Gene Scoring / Auto Cell Typing / Deconvolution | No |
| **4** | Manual Cell Typing (skippable) | Yes (iterative with `celltype_map.json`) |
| **5** | Downstream Biology | No |
| **6–7** | Meta Analysis (optional) | No |

See root `README.md` (Pipeline Workflow section) for the full Mermaid diagram and checkpoint filenames.

---

## Running a Project Pipeline

From the repo root:

```bash
snakemake -d projects/visium_hd/robin -s projects/visium_hd/robin/Snakefile all
```

Or from the project directory:

```bash
cd projects/visium_hd/robin && snakemake -d . -s Snakefile all
```

See [project_setup.md](../docs/project_setup.md) for container and environment details.

---

## Testing

Each project has `tests/` for integration and smoke tests:

```bash
pytest projects/visium/ggo_visium/tests/ -v
```

Implementation order: (1) project tests, (2) sc_tools package tests (`sc_tools/tests/`), (3) implement functions.
