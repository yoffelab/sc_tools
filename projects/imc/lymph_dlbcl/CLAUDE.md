# lymph_dlbcl IMC — Claude Code Configuration

## Sync Before Work

1. @Mission.md — current todo list and phase status
2. @Journal.md — recent decisions and conventions

For repo-wide rules (container, conventions, testing): see repo root CLAUDE.md.

---

## Project Context

**Goal:** Reproduce all figures from the DLBCL IMC manuscript (Integrative spatial analysis
reveals a hierarchy of cellular organization in diffuse large B-cell lymphoma, v7.4).
328 treatment-naive tumors, 52 markers across 2 IMC panels (immune + stromal), 12 major
cell types -> 30 subpopulations, 5 LME (lymphoma microenvironment) classes.

**Platform:** IMC (Imaging Mass Cytometry)

**Runtime:** cayuga (`/home/fs01/juk4007/elementolab/sc_tools/projects/imc/lymph_dlbcl`)

**Phase status** (see Mission.md for full detail):
- `ingest_load` / `qc_filter` / `metadata_attach`: Complete (Seurat-converted h5ad objects reused)
- `preprocess` / `celltype_manual`: Complete (existing labels from DLC_code)
- `scoring`: Complete (LME classes computed)
- **`biology` (Active):** Manuscript figure reproduction (64 PDFs generated, visual QA in progress)

---

## Running This Project

```bash
# From repo root:
./scripts/run_container.sh projects/imc/lymph_dlbcl python projects/imc/lymph_dlbcl/scripts/<script>.py
snakemake -d projects/imc/lymph_dlbcl -s projects/imc/lymph_dlbcl/Snakefile <target>

# Common targets:
snakemake -d projects/imc/lymph_dlbcl -s projects/imc/lymph_dlbcl/Snakefile fig2_lme_classes
snakemake -d projects/imc/lymph_dlbcl -s projects/imc/lymph_dlbcl/Snakefile all

# Tests:
pytest projects/imc/lymph_dlbcl/tests/ -v
```

---

## Key Files and Checkpoints

| Path | Description |
|------|-------------|
| `results/adata.celltyped.h5ad` | Primary analysis input: immune + stromal panels merged |
| `results/seurat_converted/` | Source h5ad objects converted from Seurat (47 files on cayuga) |
| `metadata/sample_metadata.csv` | 328 patient clinical metadata |
| `metadata/gene_signatures.json` | IMC marker signatures for LME class scoring |
| `figures/manuscript/` | Reproduced manuscript figures |
| `manuscript/Figures_v7.1/` | Adobe Illustrator originals for visual comparison |

---

## Project-Specific Conventions

- **Panels:** T2 (immune) + S2 (stromal) kept separate; merge only for cross-panel analyses
- **Cell types:** 12 major -> 30 subpopulations; labels from DLC_code (pre-existing)
- **LME classes:** 5 classes per manuscript; stored in `obs['LME_class']`
- **Reference data (read-only):** `/home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2`
- **Statistical correction:** Benjamini-Hochberg (FDR) for all multiple comparisons
- **Figure format:** Match manuscript figure layout exactly for visual QA comparison
