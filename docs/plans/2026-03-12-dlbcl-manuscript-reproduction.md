---
status: in_progress
created: 2026-03-12
summary: "DLBCL IMC manuscript reproduction: all 13 figures (5 main + 8 supp) from Integrative spatial analysis paper (v7.4), 328 tumors, 2 panels, 5 LME classes"
---

# DLBCL IMC Manuscript Reproduction

**Project:** `projects/imc/lymph_dlbcl`
**Runtime:** cayuga (`/athena/cayuga_0014/scratch/juk4007/sc_tools/projects/imc/lymph_dlbcl`)
**Reference:** manuscript `DLBCL_hyperion_draft_v7.4_comments.docx`, original figures `manuscript/Figures_v7.1/`

## Objective

Reproduce all 13 figures (5 main + 8 supplementary) from:
> "Integrative spatial analysis reveals a hierarchy of cellular organization in diffuse large B-cell lymphoma"

328 treatment-naive tumors, 52 markers, 2 IMC panels (immune T2 + stromal S2), 12 major cell types → 30 subpopulations, 5 LME classes.

## Data Sources

| Object | Path | Description |
|--------|------|-------------|
| Immune p4 | `results/adata.immune.celltyped.p4.h5ad` | 1.63M cells, 49 markers, 138 samples, 19 cell types. **X = zeros; real data in `layers['raw']`** |
| Stromal p4 | `results/adata.stromal.celltyped.p4.h5ad` | 1.55M cells, 50 markers. Stromal has 1 sample issue (orig.ident not mapped) |
| Spatial communities | `results/seurat_converted/stroma_spatial/SO_k30_community_cluster.h5ad` | 1.65M cells, 30 community types |
| TME object | `results/seurat_converted/s_tme_seurat.h5ad` | 576 obs (sample-level), 34 vars |
| Clinical | `metadata/DLC380_clinical.tsv` | 348 cases; FINAL_COHORT=YES → 332; columns: DLC_ID, OS_time/event, PFS_time/event, LYMPH2CX_COO |
| LME assignments | `metadata/lme_class_assignments.csv` | Sample → LME_class mapping (5 classes) |

## Key Engineering Decisions

- **X=zeros fix:** All scripts use `adata.layers.get("raw", adata.X) if adata.layers else adata.X` pattern for expression data
- **Cell type column priority:** `cluster` > `celltype` > `celltype_broad`; skip columns with all-Unknown values
- **COO filtering:** Always exclude "FAIL", "NA", "NAN" from COO values before any analysis
- **LME join:** Normalize DLC IDs via `DLC[_\s-]?(\d+)` → `DLC_{N:04d}` before any merge
- **Figure infrastructure:** All scripts import from `scripts/figure_config.py` (LME_COLORS, LME_ORDER, CELLTYPE_COLORS, COO_COLORS, build_celltype_palette)

## Figure Status

### Main Figures

| Figure | Script | Status | Notes |
|--------|--------|--------|-------|
| Fig 1: Single-cell atlas | `fig1_single_cell_atlas.py` | **PASS** | UMAP (immune+stromal), heatmap (z-scored 40 markers × 16 T/M subtypes), prevalence, fig1e COO×LME (all 5 classes, chi2 p=1.35e-03) |
| Fig 2: LME classes | `fig2_lme_classes.py` | **PASS** | Heatmap, stacked bar (30 subtypes distinct colors), violin, LME proportions |
| Fig 3: Clinical | `fig3_clinical.py` | **PASS** | OS KM p=0.0147, PFS KM p=0.0024, Cox forest (Cytotoxic highest HR), COO×LME FAIL filtered |
| Fig 4: Spatial | `fig4_spatial.py` | **PASS** | Community composition 30 types, diversity bar, enrichment bubble, ML ROC+FI |
| Fig 5: CosMx | `fig5_cosmx.py` | **PASS** | 5 panels: UMAP (241K non-B cells), dotplot (24 markers), spatial CTMA121, 6 pathways, receptor-ligand |

### Supplementary Figures

| Figure | Script | Status | Notes |
|--------|--------|--------|-------|
| Supp 1: QC panels | `supp_fig1_qc_panels.py` | **PASS** | Markers, distributions, per-sample QC |
| Supp 2: B cell | `supp_fig2_bcell.py` | **PASS** | UMAP ×4, heatmap ×4 |
| Supp 3: T cell/myeloid | `supp_fig3_tcell_myeloid.py` | **PASS** | UMAP, heatmap, proportions |
| Supp 4: Vessel | `supp_fig4_vessel.py` | **PASS** | CD31 violin, vessel density by LME |
| Supp 5: TME sensitivity | `supp_fig5_tme_sensitivity.py` | **PASS** | Cluster metrics, stability, k=5/8/10 heatmaps |
| Supp 6: Mutations | `supp_fig6_mutations.py` | **BLOCKED** | Needs `CTMA121_mut_table.csv` from cayuga |
| Supp 7: RNA-protein | `supp_fig7_rna_protein.py` | **PASS** | 7 overlapping markers, per-marker scatter + correlation bar (3 common cell types) |
| Supp 8: Extended survival | `supp_fig8_extended_survival.py` | **PASS** | Multivariate Cox (Cytotoxic highest HR), KM by COO subgroup (GCB/ABC/U) |

## Remaining Work

### Immediate

1. ~~**Fig 5 (CosMx)**~~ — DONE (2026-03-13)
2. ~~**Supp 7 (RNA-protein comparison)**~~ — DONE (2026-03-13)

### Blocked

3. **Supp 6 (Mutations)** — needs `CTMA121_mut_table.csv`; check cayuga path `/home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2/`

### Validation

4. **Numerical validation (6.2)** — verify cell counts, LME proportions, KM p-values, AUC against manuscript values
5. **End-to-end Snakemake run (6.3)** — confirm full DAG from data build → all figures

## Snakemake

```bash
# From repo root (cayuga):
snakemake -d projects/imc/lymph_dlbcl -s projects/imc/lymph_dlbcl/Snakefile all --cores 8

# Individual figures:
snakemake -d projects/imc/lymph_dlbcl -s projects/imc/lymph_dlbcl/Snakefile fig3_clinical
```
