# Container images (Apptainer/Singularity)

Pipeline runs use **Apptainer** (or Singularity) for portability. The default SIF path is:

- `containers/sc_tools.sif`

## Build the SIF

From the repository root.

**Option A (from existing Docker image):**

```bash
docker build -t sc_tools:latest .
apptainer build containers/sc_tools.sif docker-daemon://sc_tools:latest
```

**Option B (if image is in a registry):**

```bash
apptainer build containers/sc_tools.sif docker://ghcr.io/your-org/sc_tools:latest
```

Override the path with `SC_TOOLS_SIF` when running `./scripts/run_apptainer.sh` or Snakemake.

See **project_setup.md** for full usage.
