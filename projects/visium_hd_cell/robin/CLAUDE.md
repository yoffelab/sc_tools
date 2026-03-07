# robin — Claude Code Configuration

## Sync Before Work

1. @Mission.md — current todo list and phase status
2. @journal_summary.md — recent decisions

For repo-wide rules (container, conventions, testing): see repo root CLAUDE.md.

---

## Project Context

**Platform:** visium_hd_cell
**Path:** `projects/visium_hd_cell/robin`

(Describe scientific objective here.)

---

## Running This Project

```bash
# From repo root:
./scripts/run_container.sh projects/visium_hd_cell/robin python scripts/<script>.py
snakemake -d projects/visium_hd_cell/robin -s projects/visium_hd_cell/robin/Snakefile <target>

# Local conda (no container):
conda activate robin
SC_TOOLS_RUNTIME=none snakemake -d projects/visium_hd_cell/robin -s projects/visium_hd_cell/robin/Snakefile <target>

# Tests:
pytest projects/visium_hd_cell/robin/tests/ -v

# Create conda env (one-time setup, from repo root):
#   conda create -n robin python=3.11 -y && conda activate robin
#   uv pip install -e ".[deconvolution]"
```

---

## Key Files

| Path | Description |
|------|-------------|
| `results/adata.raw.h5ad` | qc_filter phase output |
| `results/adata.annotated.h5ad` | metadata_attach phase output |
| `results/adata.normalized.h5ad` | preprocess phase output |
| `results/adata.scored.h5ad` | scoring phase output (primary analysis input) |
| `results/adata.celltyped.h5ad` | celltype_manual phase output |
| `metadata/gene_signatures.json` | Gene signatures |
| `figures/manuscript/` | Publication figures |
