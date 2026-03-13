---
status: complete
created: 2026-03-12
summary: "DLBCL Fig 5: CosMx SMI analysis — 350K LME cells, 1000-plex RNA, cell type atlas, community analysis, pathway scores, receptor-ligand interactions"
---

# DLBCL Fig 5: CosMx Spatial Transcriptomics

**Project:** `projects/imc/lymph_dlbcl` (figures) + `projects/cosmx_1k/lymph_dlbcl` (data/pipeline)
**Manuscript caption:** "Fig. 5: Multicellular interactions direct the organization and function of cell communities"

## Data

| File | Path | Description |
|------|------|-------------|
| RNA object | `projects/cosmx_1k/lymph_dlbcl/data/DLBCL_cosmx_run1_CTMA.RDS` | 1.1 GB; 464,797 cells × 1,000 RNA genes; 2 tissues (CTMA121, CTMA100); clusters a–j (RNA_nbclust) |
| PRT object | `projects/cosmx_1k/lymph_dlbcl/data/DLBCL_cosmx_run1_CTMA_PRT.RDS` | 4.8 GB; 687,520 cells × 67 protein markers; has `final_cell_type` (15 types); UMAP; spatial coords (x_slide_mm, y_slide_mm) |

### PRT Cell Types (protein panel, annotated)

adipose, bcell, cd4_tcell, cd8_tcell, dendritic, endothelial, epithelial, fibroblast, immune, macrophage_monocyte, neutrophil, nk, tcell, Unknown, vascular_smooth_muscle

### Key Metadata Columns

- `final_cell_type` — cell type label (PRT only)
- `x_slide_mm`, `y_slide_mm` — spatial coordinates
- `Run_Tissue_name` — sample ID (CTMA121_PRT, CTMA100_PRT)
- `fov` — field of view

## Manuscript Figure Panels

From `Fig. 5` caption (v7.4):

| Panel | Description | Source |
|-------|-------------|--------|
| **5a** | Schematic of CosMx SMI imaging approach | Placeholder (BioRender original) |
| **5b** | UMAP of ~350K LME cells colored by cell type (B cells omitted) | RNA object (1000-plex), clusters annotated |
| **5c** | Dotplot of cell populations × key marker genes | RNA object (1000-plex) |
| **5d** | Community analysis — similar multicellular patterns as IMC | RNA object spatial clusters |
| **5e** | Pathway scores per community (TNFa, NFkb, JAK-STAT, TGFb, VEGF) | RNA object, Hallmark/pathway scoring |
| **5f** | Receptor-ligand interactions — CD8 T cell signaling across LME contexts | RNA object, CellChat or NicheNet |

## Implementation Plan

### Step 1: Convert RDS to h5ad

Script: `projects/cosmx_1k/lymph_dlbcl/scripts/convert_cosmx_to_h5ad.R`

- PRT → `results/cosmx_prt.h5ad` (protein panel, cell types, UMAP, spatial)
- RNA → `results/cosmx_rna.h5ad` (1000-plex RNA, cluster labels, spatial)

Using `SeuratDisk::SaveH5Seurat` + `Convert()`. Key things to preserve:
- PRT: `final_cell_type`, `approximateumap_*` reduction, `x_slide_mm/y_slide_mm`
- RNA: cluster labels (a–j), UMAP reduction, spatial coordinates

### Step 2: Annotate RNA clusters → cell types

Script: `projects/cosmx_1k/lymph_dlbcl/scripts/annotate_rna_clusters.py`

- Score each cluster (a–j) against lineage marker gene sets
- Map to cell types consistent with PRT labels and manuscript Supp Fig 13
- Key markers:
  - B cells: CD19, CD20, MS4A1, PAX5, BACH2
  - CD4 T: CD4, FOXP3, CXCR5, IL2RA
  - CD8 T: CD8A, CD8B, GZMB, PRF1, LAG3
  - NK: NCAM1, KLRD1, NKG7
  - Macrophage/Monocyte: CD68, CD163, MRC1, CSF1R, CD14
  - Dendritic: ITGAX, CLEC9A, CLEC4C
  - Fibroblast/CAF: PDPN, FAP, COL1A1, ACTA2
  - Endothelial: PECAM1, VWF, CDH5
  - B cells to exclude: BCL2, BCL6, MYC

### Step 3: Generate Fig 5 panels

Script: `projects/imc/lymph_dlbcl/scripts/fig5_cosmx.py`

Output: `figures/manuscript/fig5/`

#### Panel 5b — UMAP colored by cell type
- Load RNA h5ad (or PRT h5ad as fallback)
- Filter out B cells / bcell / lymphoma
- UMAP colored by annotated cell type
- Use CELLTYPE_COLORS from figure_config
- Rasterize points (PDF size)

#### Panel 5c — Dotplot (populations × marker genes)
- Key lineage marker genes from 1000-plex panel
- Mean expression + % expressing per cell type
- Z-score or arcsinh normalization
- Order cell types by lineage (T cells, B cells, myeloid, stromal, endothelial)

#### Panel 5d — Spatial map
- x_slide_mm, y_slide_mm scatter, colored by cell type
- One representative sample (CTMA121)
- Subsample if >200K cells for PDF size

#### Panel 5e — Pathway scores per community
- Score Hallmark gene sets (TNFa, NFkb, JAK-STAT, TGFb, VEGF, Hypoxia)
- Aggregate per spatial cluster / community
- Heatmap: communities × pathways (z-scored)

#### Panel 5f — Receptor-ligand (CD8 context)
- CellChat or simple ligand-receptor scoring
- Focus on CD8 T cell incoming signals across LME communities

### Step 4: Supp Fig 7 (RNA-protein comparison)

Script: `projects/imc/lymph_dlbcl/scripts/supp_fig7_rna_protein.py`

- For each overlapping cell type: correlate mean IMC protein expression (from p4 h5ad) vs mean CosMx RNA expression (from cosmx_rna.h5ad)
- 15–20 overlapping markers (CD3, CD4, CD8, CD20, CD68, CD163, PDPN, CD31, Ki67, BCL2, BCL6, etc.)
- Scatter plot per marker: IMC protein intensity vs CosMx RNA counts (Pearson r)
- Summary correlation heatmap across all markers

## File Outputs

```
projects/cosmx_1k/lymph_dlbcl/
  scripts/
    convert_cosmx_to_h5ad.R        # RDS → h5ad conversion (Step 1)
    annotate_rna_clusters.py       # cluster → cell type annotation (Step 2)
  results/
    cosmx_prt.h5ad                 # protein panel h5ad
    cosmx_rna.h5ad                 # RNA panel h5ad
    cosmx_rna_annotated.h5ad       # RNA panel with cell type labels

projects/imc/lymph_dlbcl/
  scripts/
    fig5_cosmx.py                  # main figure script (Step 3)
    supp_fig7_rna_protein.py       # RNA-protein correlation (Step 4)
  figures/manuscript/
    fig5/                          # Fig 5 panels
    supp_fig7/                     # Supp Fig 7 panels
```

## Blockers / Dependencies

- [x] CosMx data available locally (both RDS files present)
- [ ] SeuratDisk conversion — test on PRT first (~10 min for 4.8 GB RDS)
- [ ] RNA cluster annotation — requires scoring against known markers
- [ ] Pathway scoring — requires Hallmark gene sets (available via `sc_tools.tl.load_hallmark()`)
- [ ] CellChat — requires R or Python cellchat equivalent; may simplify to NicheNet scores
- [ ] Supp Fig 7 — depends on CosMx RNA h5ad + IMC p4 h5ad

## Success Criteria

1. Fig 5b UMAP shows clearly separated non-B cell populations (same cell types as PRT labels)
2. Fig 5c dotplot shows expected lineage marker patterns (CD8A in CD8 T cells, PDPN in fibroblasts, etc.)
3. Fig 5d spatial map shows spatially coherent cell type distributions
4. Panel count matches manuscript (at minimum 5b + 5c produced; 5d-5f best effort)
5. Cell count in UMAP legend ~350K non-B cells (matching manuscript)
