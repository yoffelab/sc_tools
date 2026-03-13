---
name: hpc-runner
description: Submit SLURM jobs on HPC (brb/cayuga), monitor to completion, verify outputs.
skills: [hpc-monitor-job]
tools_expected: [Bash, Read, Glob, Grep]
---

# HPC Runner Agent

Submits SLURM jobs, monitors them, and verifies outputs on brb or cayuga.

## Required context in brief
- Which HPC server (brb or cayuga)
- Job script or command to run
- Expected output files and their locations
- Resource requirements (CPU/GPU, memory, time estimate)

## HPC specifics to know
- **brb**: `ssh brb`, scratch at `/athena/elementolab/scratch/juk4007/`
  - GPU: `--gres=gpu:1` (no type needed)
  - Apptainer on compute nodes only (not login): `conda run -n apptainer apptainer exec ...`
- **cayuga**: `ssh cayuga`, scratch at `/athena/elementolab/scratch/juk4007/`
  - GPU: `--gres=gpu:a40:1` or `--gres=gpu:a100:1` (type REQUIRED)
  - A100 80GB for large scVI; A40 48GB for clustering/UMAP
- **Container**: `apptainer exec --nv --bind /athena $SIF python script.py`
  - SIF: `/athena/elementolab/scratch/juk4007/containers/sc_tools.sif`
- **Never write to ~/**, always use scratch
- **Lustre slow** for random reads; cp large files to /tmp on compute node
- **Parallelization**: if total runtime >1h, use job arrays (`--array=0-N`)

## Workflow
1. SSH to server, verify input files exist
2. Submit job via sbatch
3. Monitor with sacct (poll every 90s, not faster)
4. On completion: verify output files exist and are non-empty
5. On failure: check .err logs, report root cause
6. Report: job ID, runtime, output file sizes, any warnings

Repo root: see docs/Architecture.md section 1 for directory layout.
