# sc_tools — Spatial and Single-Cell Omics Toolkit

[![CI](https://github.com/yoffelab/sc_tools/actions/workflows/ci.yml/badge.svg)](https://github.com/yoffelab/sc_tools/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/sci-sc-tools)](https://pypi.org/project/sci-sc-tools/)
[![Python](https://img.shields.io/pypi/pyversions/sci-sc-tools)](https://pypi.org/project/sci-sc-tools/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Docs](https://readthedocs.org/projects/sc-tools/badge/?version=latest)](https://sc-tools.readthedocs.io)

**sc_tools** is a reusable Python toolkit and pipeline for spatial transcriptomics and single-cell omics. It wraps scanpy, squidpy, scVI-tools, and related libraries in a consistent API and provides a phased, Snakemake-driven workflow for multi-modality projects (Visium, Visium HD, Xenium, IMC, CosMx).

The package is installable via `pyproject.toml` and designed to be shared across projects. Project-specific data, results, and figures all live under `projects/<platform>/<project_name>/`.

---

## Installation

From the repository root:

```bash
# Minimal install (plotting, QC, gene scoring)
pip install -e .

# Full pipeline (deconvolution: scvi-tools, tangram)
pip install -e ".[deconvolution]"

# Common combinations
pip install -e ".[deconvolution,integration,geneset]"
pip install -e ".[deconvolution,gpu]"
```

**Available extras:**

| Extra | Installs | When to use |
|-------|----------|-------------|
| `deconvolution` | scvi-tools, tangram-sc | Cell-type deconvolution (Cell2location, Tangram, DestVI) |
| `integration` | harmonypy | Harmony batch correction |
| `geneset` | gseapy, pyucell | GSEA pseudobulk, UCell scoring |
| `decoupler` | decoupler | TF/pathway activity |
| `spatial` | utag | Spatial-aware clustering (UTAG) |
| `gpu` | torch, rapids-singlecell | GPU-accelerated preprocessing |
| `viz` | marsilea | Declarative composite figures |
| `benchmark` | scikit-image, scib-metrics, cellpose, stardist | Segmentation/integration benchmarking |
| `dev` | pytest, ruff | Development and testing |
| `docs` | sphinx, pydata-sphinx-theme, myst-nb | Build documentation |

**Lint (required before commit):** `make lint` runs Ruff check and format on `sc_tools`.

**Container:** See [project_setup.md](project_setup.md) for Apptainer/Docker setup and per-project usage.

---

## What sc_tools provides

| Module | Purpose |
|--------|---------|
| `sc_tools.pl` | Spatial plots, heatmaps, statistical annotations, volcano plots, GSEA dotplots, versioned PDF+PNG saving |
| `sc_tools.tl` | Signature scoring (scanpy/UCell/ssGSEA), gene set loaders (Hallmark bundled), ORA/GSEA, colocalization, deconvolution |
| `sc_tools.qc` | Per-sample QC filtering, cross-sample comparison plots, HTML QC reports |
| `sc_tools.pp` | Modality-aware preprocessing recipes: normalize, integrate (scVI/Harmony/CytoVI), reduce, cluster |
| `sc_tools.ingest` | Batch manifests (`metadata/phase0/`), Space Ranger/Xenium/IMC command builders, AnnData loaders per modality, `concat_samples()` |
| `sc_tools.validate` | Checkpoint validation for p1–p4 (required obs keys, obsm, representation); auto-fix; CLI for Snakemake |
| `sc_tools.memory` | Memory profiling, GPU detection and auto-backend selection |

---

## Pipeline Workflow

The pipeline is **non-linear** with human-in-loop phases. Branching points and explicit input files (e.g. clinical metadata) bypass manual steps.

```mermaid
flowchart LR
    subgraph P0["Phase 0: Upstream Data Processing"]
        direction TB
        Z1["(0a) HPC: Space Ranger / Xenium Ranger / IMC pipeline"] --> Z2["(0b) Load per sample → adata.p0.h5ad"]
    end

    subgraph P1["Phase 1: QC and Concatenation"]
        direction TB
        A1[Load Phase 0 AnnData] --> A2[Per-sample QC + filter]
        A2 --> A3[concat_samples()]
        A3 --> A4["QC report: figures/QC/raw/"]
    end

    subgraph P2["Phase 2: Metadata Attachment"]
        direction TB
        B1{Clinical map provided?} --> |"No: HIL"| B2[Prepare CSV/xlsx]
        B2 --> B1
        B1 --> |"Yes"| B3[Join to adata.obs]
    end

    subgraph P3["Phase 3: Preprocessing"]
        direction TB
        C1[Backup adata.raw] --> C2[Filter QC failures]
        C2 --> C3[Normalize, batch correct, cluster]
        C3 --> C5["QC report: figures/QC/post/"]
    end

    subgraph P35["Phase 3.5: Demographics"]
        direction TB
        D1[Cohort stats]
        D2[Figure 1]
    end

    subgraph P35b["Phase 3.5b: Gene Scoring, Cell Typing, Deconvolution"]
        direction TB
        D3[Hallmark + project sigs in obsm]
        D4[Automated cell typing]
        D5[Cell-type deconvolution optional]
    end

    subgraph P4["Phase 4: Manual Cell Typing"]
        direction TB
        E1[Extract cluster_id] --> E2[JSON: cluster_id→celltype]
        E2 --> E3{Satisfactory?}
        E3 --> |"No: HIL"| E2
        E3 --> |"Yes"| E4[Apply celltype + celltype_broad]
    end

    subgraph P5P67["Phase 5 & 6–7"]
        direction TB
        subgraph P5["Phase 5: Downstream Biology"]
            direction TB
            F1[Spatial/process analysis]
            F2[Colocalization, Moran's I]
            F3[Publication figures]
        end
        subgraph P67["Phase 6–7: Meta Analysis"]
            direction TB
            G1[Phase 6: Aggregate ROI/patient]
            G2[Phase 7: Downstream on aggregated]
        end
    end

    P0 --> P1 --> P2 --> P3
    P3 --> P35
    P3 --> P35b
    P35b --> P4 --> P5P67
    P35b -.-> |"Skip Phase 4 if automated typing adequate"| P5P67
    P1 -.-> |"START HERE: processed outs available"| P2
    P2 -.-> |"START HERE: preprocessed AnnData available"| P3
    P3 -.-> |"START HERE: clustered AnnData available"| P35b
    P35b -.-> |"START HERE: phenotyped AnnData available"| P4

    style P0 fill:#f0f4ff
    style P1 fill:#e3f2fd
    style P2 fill:#fff3e0
    style P3 fill:#e8f5e9
    style P35 fill:#e1f5fe
    style P35b fill:#e8eaf6
    style P4 fill:#fff3e0
    style P5 fill:#f3e5f5
    style P67 fill:#fafafa
```

| Phase | Name | Human-in-Loop? | Checkpoint |
|-------|------|----------------|------------|
| **0a** | Platform tools (Space Ranger / Xenium / IMC) | No | `data/{sample_id}/outs/` or `processed/{sample}/` |
| **0b** | Load per-sample into AnnData | No | `data/{sample_id}/adata.p0.h5ad` |
| **1** | QC and Concatenation | No | `results/adata.raw.p1.h5ad` |
| **2** | Metadata Attachment | Yes (unless map provided) | `results/adata.annotated.p2.h5ad` |
| **3** | Preprocessing | No | `results/adata.normalized.p3.h5ad` |
| **3.5** | Demographics (branch) | Project-specific | Figure 1 |
| **3.5b** | Gene Scoring / Auto Cell Typing / Deconvolution | No | `results/adata.normalized.scored.p35.h5ad` |
| **4** | Manual Cell Typing (skippable) | Yes (iterative) | `results/adata.celltyped.p4.h5ad` |
| **5** | Downstream Biology | No | `figures/manuscript/` |
| **6–7** | Meta Analysis (optional) | No | `results/adata.{level}.{feature}.h5ad` |

---

## Repository layout

```text
sc_tools/               # Reusable Python package (pl, tl, qc, pp, ingest, validate, memory, utils)
scripts/                # Shared scripts and run_container.sh
projects/               # All project-specific content
    create_project.sh   # ./projects/create_project.sh <project_name> <data_type>
    <platform>/
        <project_name>/
            data/           # Raw sequencing and imaging
            metadata/       # Gene signatures, sample_metadata.csv, celltype_map.json
            results/        # AnnData checkpoints (.h5ad)
            figures/        # QC, manuscript figures
            scripts/        # Project-specific analysis scripts
            outputs/        # Intermediate outputs (deconvolution logs, etc.)
            tests/          # Project integration tests
            Snakefile       # Project pipeline
            Mission.md      # Project todo list
            Journal.md      # Project decision log
docs/                   # Sphinx API documentation (make docs)
containers/             # Apptainer SIF and Dockerfile
```

Key docs: **`Mission.md`** (toolkit todos), **`Architecture.md`** (data flow, checkpoint contracts), **`skills.md`** (statistical and coding standards), **`project_setup.md`** (container and environment setup).

---

## Running the pipeline

Each project uses **Snakemake** as the workflow engine. From the repo root:

```bash
# Run a project target
snakemake -d projects/visium_hd/robin -s projects/visium_hd/robin/Snakefile all

# Run with container (auto-detects Apptainer on Linux, Docker on macOS)
./scripts/run_container.sh projects/visium_hd/robin python scripts/run_signature_scoring.py

# Override runtime
SC_TOOLS_RUNTIME=docker ./scripts/run_container.sh projects/visium/ggo_visium
```

See [project_setup.md](project_setup.md) for container build instructions and per-project setup.

---

## License and contact

Internal research project. For questions or collaboration, contact the repository maintainers or the Yoffe Lab.
