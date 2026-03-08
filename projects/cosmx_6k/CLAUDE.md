# cosmx_6k — Claude Code Configuration

## Sync Before Work

1. @Mission.md — current todo list and phase status
2. @journal_summary.md — recent decisions

For repo-wide rules (container, conventions, testing): see repo root CLAUDE.md.

---

## Project Context

**Platform:** cosmx
**Path:** `projects/cosmx_6k`

(Describe scientific objective here.)

---

## Running This Project

```bash
# From repo root:
./scripts/run_container.sh projects/cosmx_6k python scripts/<script>.py
snakemake -d projects/cosmx_6k -s projects/cosmx_6k/Snakefile <target>

# Local conda (no container):
conda activate cosmx_6k
SC_TOOLS_RUNTIME=none snakemake -d projects/cosmx_6k -s projects/cosmx_6k/Snakefile <target>

# Tests:
pytest projects/cosmx_6k/tests/ -v

# Create conda env (one-time setup, from repo root):
#   conda create -n cosmx_6k python=3.11 -y && conda activate cosmx_6k
#   uv pip install -e ".[deconvolution]"
```

---

## Key Files

| Path | Description |
|------|-------------|
| `results/adata.filtered.h5ad` | qc_filter output |
| `results/adata.scored.h5ad` | scoring output (primary analysis input) |
| `metadata/gene_signatures.json` | Gene signatures |
| `figures/manuscript/` | Publication figures |
