---
status: in_progress
created: 2026-03-12
summary: Build and deploy sc_tools Apptainer SIF to HPC servers (brb + cayuga) to eliminate per-job conda dependency installation pain.
---

# Container Deployment Plan (v1)

## Goal

Replace ad-hoc `conda activate sc_tools` + manual `pip install` in SLURM jobs with a frozen
Apptainer SIF that contains all pipeline dependencies. Build once locally → deploy to both HPC
servers → all future jobs use `apptainer exec --nv $SIF python script.py`.

## Why

Repeated SSH sessions to install packages (`lifelines`, `leidenalg`, `scanorama`, etc.) were the
main source of HPC friction. Conda envs on brb and cayuga drift independently and are painful to
sync. The container fixes this with a single frozen image.

## What Was Added

### `pyproject.toml` — new `[pipeline]` extra

Bundles everything needed for HPC runtime into one install target:
- `scvi-tools`, `tangram-sc` (deconvolution)
- `harmonypy`, `bbknn`, `scanorama` (integration benchmarks)
- `leidenalg`, `igraph` (leiden clustering — was buried in `dev`)
- `lifelines>=0.27` (survival analysis — was missing entirely, caused cayuga pain)
- `scib-metrics`, `plotly` (integration benchmarking reports)
- `sqlalchemy`, `alembic` (registry)
- `marsilea`, `fsspec`, `s3fs`, `zarr` (viz + storage)

### `Dockerfile`

- Installs `gfortran` (needed for `scikit-misc` from `scanpy[skmisc]`)
- Installs CUDA 12.1 torch first (so scvi uses GPU when `--nv` is passed)
- Installs `.[pipeline]` extra
- Verifies `sc_tools`, `scvi`, `lifelines`, `harmonypy` on build

### `Makefile` — container targets

```makefile
make container-build          # docker buildx build --platform linux/amd64 --output type=docker
make container-test           # local smoke test
make container-sif-via-brb   # OCI export → rsync to brb → apptainer build → rsync SIF back
make container-deploy         # rsync SIF to brb + cayuga, verify imports
make container-all-macos      # full macOS workflow: build + sif-via-brb + deploy
```

### `containers/sc_tools.sif` (gitignored, 3.7GB)

Apptainer SIF deployed to:
- `brb:/athena/elementolab/scratch/juk4007/containers/sc_tools.sif` ✅ built and verified
- `cayuga:/athena/elementolab/scratch/juk4007/containers/sc_tools.sif` — transfer in progress

## Current Status (2026-03-12)

| Step | Status | Notes |
|------|--------|-------|
| `pyproject.toml` — `[pipeline]` extra | ✅ Done | |
| `Dockerfile` — gfortran + CUDA torch | ✅ Done | |
| `Makefile` — container targets | ✅ Done | |
| Docker image build (`linux/amd64`, 4GB) | ✅ Done | `sc_tools:latest` |
| OCI tar export (3.8GB) | ✅ Done | `/tmp/sc_tools_oci.tar` |
| Apptainer installed on brb via conda | ✅ Done | `conda env apptainer`, apptainer 1.4.5 |
| SIF built on brb (SLURM job 13991073) | ✅ Done | 6m58s, 3.7GB |
| SIF verified on brb compute (job 13992017) | ✅ Done | `brb compute: all imports OK` |
| SIF rsynced to local `containers/` | ✅ Done | 3.7GB |
| SIF deployed to cayuga | 🔄 In progress | rsync running now |
| SIF verified on cayuga | ⬜ Pending | |
| SLURM scripts updated to use container | ⬜ Pending | See next steps below |

## Limitation: brb Login Node

`apptainer exec` fails on the brb **login node** (`user namespace disabled`, CentOS 7 kernel).
Works fine on compute nodes (confirmed in job 13992017). Use `apptainer` only inside SLURM jobs
on brb; the login node is login-only anyway.

## Next Steps

### Immediate (after cayuga verify)

1. **Verify SIF on cayuga** — run a quick test job:
   ```bash
   ssh cayuga "module load apptainer && apptainer exec \
       /athena/elementolab/scratch/juk4007/containers/sc_tools.sif \
       python -c 'import sc_tools, scvi, lifelines; print(\"cayuga: OK\")'"
   ```

2. **Update existing SLURM scripts** to use the container. Priority order:
   - `projects/multiplatform/ibd_spatial/scripts/run_m1.sh` (cayuga GPU, scVI)
   - `projects/multiplatform/ibd_spatial/scripts/run_m2.sh` (cayuga GPU, scVI)
   - `projects/visium_hd_cell/robin/scripts/run_pipeline.sbatch` (brb GPU, scVI)
   - `projects/imc/lymph_dlbcl/scripts/run_on_cayuga.sh`

   **Pattern change** (per script):
   ```bash
   # Before:
   eval "$(conda shell.bash hook 2>/dev/null)"
   conda activate sc_tools

   # After:
   SIF=/athena/elementolab/scratch/juk4007/containers/sc_tools.sif
   # Then prefix every python call:
   apptainer exec --nv --bind /athena "$SIF" python scripts/foo.py
   ```

3. **Update `create_project.sh`** SLURM template to default to apptainer pattern

### Future (when updating dependencies)

```bash
# On macOS, from repo root:
# 1. Edit pyproject.toml to add new package
# 2. Rebuild and redeploy:
make container-all-macos
```

This is the only step needed — no SSH, no conda install, no env drift.

### Known Issues / Watch Points

- `run_container.sh` already handles apptainer auto-detection (already correct)
- The `apptainer` conda env on brb (for building) is separate from `sc_tools` conda env
- SIF bind path: always `--bind /athena` so project data dirs are accessible
- GPU: always use `--nv` flag for scVI jobs; harmless on CPU nodes
- The local `containers/sc_tools.sif` (3.7GB) is gitignored — never commit it
