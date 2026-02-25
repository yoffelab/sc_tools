# Project Setup: Apptainer + Snakemake + Conda

This document describes how to run the sc_tools pipeline with **Apptainer** (default), build the container image, and use Snakemake or the conda environment for each project.

---

## 1. Overview

- **Default runtime: Apptainer** (Singularity-compatible) for portability (HPC, no Docker daemon). One SIF image for all projects.
- **Snakemake** is the canonical pipeline; each project has a `Snakefile` and `config.yaml`. Makefiles are deprecated (see deprecation note in each Makefile).
- **Image definition:** The repo `Dockerfile` defines the image (conda env `sc_tools`, sc_tools installed via UV). Build a Docker image first, then convert to SIF for Apptainer, or pull from a registry.
- **One image for all projects;** the project is chosen via `run_apptainer.sh <project_path>` or Snakemake `-d <project_dir>`.

---

## 2. Prerequisites

- **Apptainer** (or Singularity) installed
  - Linux: `sudo apt install apptainer` (Ubuntu) or see https://apptainer.org/docs/admin/latest/installation.html
  - macOS: Apptainer is Linux-focused; use a VM or HPC for SIF runs, or use local conda (see Local Development below)
- **Optional:** Docker (to build the image and convert to SIF via `docker-daemon://`)
- **Optional:** Conda or miniconda for local development without containers

---

## 3. Build the SIF (Apptainer image)

The SIF is built from the Dockerfile. Default path: `containers/sc_tools.sif`.

**Option A — From existing Docker image (recommended):**

From the repository root:

```bash
docker build -t sc_tools:latest .
apptainer build containers/sc_tools.sif docker-daemon://sc_tools:latest
```

**Option B — From a registry:**

If the image is pushed to a registry (e.g. GHCR):

```bash
apptainer build containers/sc_tools.sif docker://ghcr.io/your-org/sc_tools:latest
```

Override the SIF path with the environment variable `SC_TOOLS_SIF` when running scripts or Snakemake.

See **containers/README.md** for more details.

---

## 4. Running with Apptainer

### One-off script (run_apptainer.sh)

From the repository root:

```bash
./scripts/run_apptainer.sh <project_path> [command]
```

Examples:

- Interactive bash in project dir:
  ```bash
  ./scripts/run_apptainer.sh projects/visium/ggo_visium
  ```

- Run a script:
  ```bash
  ./scripts/run_apptainer.sh projects/visium/ggo_visium python scripts/score_gene_signatures.py
  ```

- Override SIF path:
  ```bash
  SC_TOOLS_SIF=/path/to/sc_tools.sif ./scripts/run_apptainer.sh projects/visium/ggo_visium python scripts/foo.py
  ```

| Project    | Path                         | Example |
|------------|------------------------------|---------|
| ggo_visium | projects/visium/ggo_visium   | `./scripts/run_apptainer.sh projects/visium/ggo_visium python scripts/score_gene_signatures.py` |
| robin      | projects/visium_hd/robin     | `./scripts/run_apptainer.sh projects/visium_hd/robin python scripts/run_signature_scoring.py` |
| ggo-imc    | projects/imc/ggo-imc         | `./scripts/run_apptainer.sh projects/imc/ggo-imc python scripts/download_yaml.py` |

### Snakemake (canonical pipeline)

From the repository root, run a project pipeline:

```bash
snakemake -d projects/visium/ggo_visium -s projects/visium/ggo_visium/Snakefile all
snakemake -d projects/visium/ggo_visium -s projects/visium/ggo_visium/Snakefile localization_only
```

Or from the project directory:

```bash
cd projects/visium/ggo_visium
snakemake -d . -s Snakefile all
```

Snakemake runs each rule by invoking `run_apptainer.sh` (as defined in each project Snakefile). The project `config.yaml` sets `repo_root` and `project_rel` so the runner finds the repo and SIF.

**Repo-root targets (lint, format, test):**

```bash
snakemake -s Snakefile lint
snakemake -s Snakefile test
snakemake -s Snakefile all
```

---

## 5. Per-Project Paths

`run_apptainer.sh` binds the full repository to `/workspace` and sets the working directory to the project (e.g. `/workspace/projects/visium/ggo_visium`). Run from the repository root so paths resolve correctly. For ggo_visium, `GGO_VISIUM_PROJECT_DIR` is set inside the container so scripts resolve project paths.

---

## 6. Local Development (Non-Container)

### Conda

```bash
conda env create -f environment.yml
conda activate sc_tools
pip install -e ".[deconvolution]"
```

Or create manually:

```bash
conda create -n sc_tools python=3.10 -y
conda activate sc_tools
pip install -e ".[deconvolution]"
```

### UV (without conda)

```bash
uv venv
source .venv/bin/activate
uv pip install -e ".[deconvolution]"
```

For Snakemake without Apptainer, you would run rules locally (edit the Snakefile to call `python scripts/...` instead of `run_apptainer.sh` when not using the container).

---

## 7. Docker (legacy / build only)

Docker is still used to **build** the image that is converted to SIF. For running pipelines, use Apptainer and Snakemake.

If you prefer to run with Docker (e.g. on macOS where Apptainer is not native):

```bash
./scripts/run_docker.sh <project_path> [command]
```

See the script header and in-repo comments. The canonical workflow is Apptainer + Snakemake.

---

## 8. Troubleshooting

**SIF not found**

Ensure the SIF is built and path is correct: `containers/sc_tools.sif` from repo root, or set `SC_TOOLS_SIF`.

**Dask: "The legacy implementation is no longer supported"**

Set before running:

```bash
export DASK_DATAFRAME__QUERY_PLANNING=true
```

**"Killed" / Exit 137 (Out of memory)**

Use a machine with more RAM or increase resource limits for the container runtime.

---

## 9. Verification

After building the SIF:

```bash
apptainer exec containers/sc_tools.sif python -c "import sc_tools; import scanpy; import squidpy; print('sc_tools OK')"
```

With the runner:

```bash
./scripts/run_apptainer.sh projects/visium/ggo_visium python -c "import sc_tools; print(sc_tools.__version__)"
```
