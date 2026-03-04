# cosmx_1k — Claude Code Configuration

## Sync Before Work

1. @Mission.md — current todo list and phase status
2. @journal_summary.md — recent decisions

For repo-wide rules (container, conventions, testing): see repo root CLAUDE.md.

---

## Project Context

**Platform:** cosmx
**Path:** `projects/cosmx/cosmx_1k`

(Describe scientific objective here.)

---

## Running This Project

```bash
# From repo root:
./scripts/run_container.sh projects/cosmx/cosmx_1k python scripts/<script>.py
snakemake -d projects/cosmx/cosmx_1k -s projects/cosmx/cosmx_1k/Snakefile <target>

# Local conda (no container):
conda activate cosmx_1k
SC_TOOLS_RUNTIME=none snakemake -d projects/cosmx/cosmx_1k -s projects/cosmx/cosmx_1k/Snakefile <target>

# Tests:
pytest projects/cosmx/cosmx_1k/tests/ -v

# Create conda env (one-time setup, from repo root):
#   conda create -n cosmx_1k python=3.10 -y && conda activate cosmx_1k
#   uv pip install -e ".[deconvolution]"
```

---

## Key Files

| Path | Description |
|------|-------------|
| `results/adata.raw.p1.h5ad` | Phase 1 output |
| `results/adata.normalized.scored.p35.h5ad` | Phase 3.5b output (primary analysis input) |
| `metadata/gene_signatures.json` | Gene signatures |
| `figures/manuscript/` | Publication figures |
