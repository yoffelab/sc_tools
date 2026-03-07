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
flowchart TD
    subgraph ING["Ingestion"]
        Z1["ingest_raw<br/>HPC: SpaceRanger / IMC"] --> CP0a[("data/{id}/outs/")]
        CP0a --> Z2["ingest_load<br/>Load per-sample → adata.h5ad"]
        Z2 --> CP0b[("adata.p0.h5ad<br/>obs: sample, library_id,<br/>raw_data_dir<br/>obsm: spatial<br/>X: raw counts")]
    end

    subgraph QC["QC and Metadata"]
        A1["qc_filter<br/>QC + Concatenation"] --> CP1[("adata.raw.h5ad<br/>obs: sample, raw_data_dir<br/>obsm: spatial<br/>X: raw counts, concatenated")]
        CP1 --> QR1[/"pre_filter_qc.html"/]
        CP1 --> A2["metadata_attach<br/>Attach clinical metadata"]
        A2 --> CP2[("adata.annotated.h5ad<br/>+ clinical metadata in obs")]
        CP2 --> QR2[/"post_filter_qc.html"/]
    end

    subgraph PRE["Preprocessing"]
        C1["preprocess<br/>Normalize + Integrate + Cluster"] --> CP3[("adata.normalized.h5ad<br/>obsm: X_scvi or embedding<br/>obs: leiden<br/>adata.raw: backed up")]
        CP3 --> QR3[/"post_integration_qc.html<br/>Batch score = primary"/]
    end

    subgraph DEM["Demographics (branch)"]
        D1["demographics<br/>Cohort stats, Figure 1"]
    end

    subgraph SCO["Scoring (branch)"]
        D3["scoring<br/>Gene scoring + Auto cell typing"] --> CP35[("adata.scored.h5ad<br/>obsm: signature_score,<br/>signature_score_z<br/>uns: signature_score_report")]
    end

    subgraph CT["Cell Typing (optional, iterative)"]
        E1["celltype_manual<br/>Cluster to celltype JSON"] --> E2{Satisfactory?}
        E2 -->|No| E1
        E2 -->|Yes| E3[Apply celltype labels]
        E3 --> CP4[("adata.celltyped.h5ad<br/>obs: celltype, celltype_broad")]
        CP4 --> QR4[/"post_celltyping_qc.html<br/>Full bio + batch scores"/]
    end

    subgraph BIO["Biology"]
        F1["biology<br/>Spatial analysis, figures"]
    end

    subgraph META["Meta Analysis (optional)"]
        G1["meta_analysis<br/>ROI / patient aggregation"] --> CP67[("adata.{level}.{feature}.h5ad<br/>obs indexed by roi or patient")]
    end

    ING --> QC --> PRE
    PRE --> DEM
    PRE --> SCO
    SCO --> CT
    SCO -.->|skip manual CT| BIO
    CT --> BIO
    BIO --> META

    START_P3(["Entry: preprocessed AnnData"]) -.-> PRE
    START_P35(["Entry: scored AnnData"]) -.-> SCO
    START_P4(["Entry: celltyped AnnData"]) -.-> CT

    style CP0a fill:#fff3e0,stroke:#ff9800
    style CP0b fill:#fff3e0,stroke:#ff9800
    style CP1 fill:#fff3e0,stroke:#ff9800
    style CP2 fill:#fff3e0,stroke:#ff9800
    style CP3 fill:#fff3e0,stroke:#ff9800
    style CP35 fill:#fff3e0,stroke:#ff9800
    style CP4 fill:#fff3e0,stroke:#ff9800
    style CP67 fill:#fff3e0,stroke:#ff9800
    style QR1 fill:#e8f4e8,stroke:#4caf50
    style QR2 fill:#e8f4e8,stroke:#4caf50
    style QR3 fill:#e8f4e8,stroke:#4caf50
    style QR4 fill:#e8f4e8,stroke:#4caf50
```

### Phase summary

| Slug | Name | Checkpoint | Required Data | QC Report |
|------|------|------------|---------------|-----------|
| `ingest_raw` | Platform tools (Space Ranger / Xenium / IMC) | `data/{sample_id}/outs/` | Platform-specific raw output | |
| `ingest_load` | Load per-sample into AnnData | `data/{sample_id}/adata.p0.h5ad` | `obs[sample, library_id, raw_data_dir]`, `obsm[spatial]`, `X` raw counts | |
| `qc_filter` | QC and Concatenation | `results/adata.raw.h5ad` | `obs[sample, raw_data_dir]`, `obsm[spatial]`, `X` raw counts, all samples concatenated | `pre_filter_qc.html` |
| `metadata_attach` | Metadata Attachment (HIL) | `results/adata.annotated.h5ad` | All of `qc_filter` + clinical columns in `obs` | `post_filter_qc.html` |
| `preprocess` | Preprocessing + Integration | `results/adata.normalized.h5ad` | `obsm[X_scvi]` (or embedding), `obs[leiden]`, `adata.raw` backed up | `post_integration_qc.html` |
| `demographics` | Demographics (branch, optional) | Figure 1 | Cohort metadata from `preprocess` | |
| `scoring` | Gene Scoring / Auto Cell Typing | `results/adata.scored.h5ad` | `obsm[signature_score, signature_score_z]`, `uns[signature_score_report]` | |
| `celltype_manual` | Manual Cell Typing (optional) | `results/adata.celltyped.h5ad` | All of `scoring` + `obs[celltype, celltype_broad]` | `post_celltyping_qc.html` |
| `biology` | Downstream Biology | `figures/manuscript/` | Reads from `scoring` or `celltype_manual` checkpoint | |
| `meta_analysis` | Meta Analysis (optional) | `results/adata.{level}.{feature}.h5ad` | `obs` indexed by roi/patient; `X` = aggregated feature | |

> Checkpoints are orange circles in the diagram; QC reports are green parallelograms. All QC reports are date-versioned (`YYYYMMDD`) under `figures/QC/`. See [Architecture.md Section 2.2](Architecture.md) for full validation contracts.

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
