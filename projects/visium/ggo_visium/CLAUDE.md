# GGO Visium — Claude Code Configuration

## Sync Before Work

1. @Mission.md — current todo list and phase status
2. @journal_summary.md — recent decisions and conventions

For repo-wide rules (container, conventions, testing): see repo root CLAUDE.md.

---

## Project Context

**Goal:** Identify transcriptional and spatial differences across lung tumor progression
(Normal → Non-Solid/GGO → Solid) with focus on TLS (Tertiary Lymphoid Structure) heterogeneity.

**Platform:** 10x Visium spatial transcriptomics

**Phase status** (see Mission.md for full detail):
- Phases 1–4: Complete
- Phase 3.5 (Demographics): Partial
- **Phase 5 (Active):** tumor_differences, macrophage localization, neutrophil–cytotoxic T-cell colocalization, process colocalization, TLS, signature heatmaps, manuscript figures
- Phases 6–7: Pending

---

## Running This Project

```bash
# From repo root:
./scripts/run_container.sh projects/visium/ggo_visium python scripts/<script>.py
snakemake -d projects/visium/ggo_visium -s projects/visium/ggo_visium/Snakefile <target>

# Common targets:
snakemake -d projects/visium/ggo_visium -s projects/visium/ggo_visium/Snakefile phase35b
snakemake -d projects/visium/ggo_visium -s projects/visium/ggo_visium/Snakefile localization_only
snakemake -d projects/visium/ggo_visium -s projects/visium/ggo_visium/Snakefile signature_heatmaps_versioned

# Tests:
pytest projects/visium/ggo_visium/tests/ -v
```

---

## Key Files and Checkpoints

| Path | Description |
|------|-------------|
| `results/adata.annotated.h5ad` | metadata_attach output: AnnData with clinical metadata |
| `results/adata.scored.h5ad` | **Primary biology input**: scores in `obsm['signature_score']` |
| `results/adata.celltyped.h5ad` | celltype_manual output: manual cell type labels |
| `results/adata.deconvolution.h5ad` | Cell-type proportions (optional) |
| `metadata/gene_signatures.json` | Gene signatures for scoring |
| `metadata/celltype_map.json` | cluster_id → celltype mapping |
| `figures/manuscript/` | Publication figures |
| `figures/QC/` | QC reports (raw/ and post/) |

---

## Project-Specific Conventions

- **Tumor stages:** Normal, Non-Solid (GGO), Solid — in that order in figures
- **Comparison mode:** 1-vs-rest default for tumor stage comparisons
- **scVI:** 30 latent dimensions
- **Deconvolution:** batch per `library_id`; process in backed mode to avoid segfaults; proportions in `obsm['cell_type_proportions']` and `obs['prop_*']`
- **Signature heatmaps:** 3-level hierarchy; solidity order Normal → Non-Solid → Solid
- **Statistical output:** significance text box in top-right corner; dodged bars for multi-group

---

## Active Scripts (Phase 5)

| Script | Output |
|--------|--------|
| `scripts/tumor_differences.py` | `figures/manuscript/tumor_differences/` |
| `scripts/macrophage_localization.py` | `figures/manuscript/macrophage_localization/` |
| `scripts/neutrophil_cytotoxic_tcell_localization.py` | `figures/manuscript/neutrophil_cytotoxic_tcell_localization/` |
| `scripts/process_colocalization.py` | `figures/process_colocalization/`, `results/process_colocalization/` |
| `scripts/signature_heatmap_versioned.py` | `figures/manuscript/signature_heatmaps_versioned.done` |
| `scripts/tls_analysis.py` | `figures/manuscript/tls_analysis/` |
