# Mission: Two-Panel IMC DLBCL (lymph_dlbcl)

**Project:** `projects/imc/lymph_dlbcl`
**Current Status:** Manuscript reproduction — building reproducible Snakemake pipeline
**Author:** Junbum Kim
**Last Updated:** 2026-03-04

This file holds **project-specific** goals. Repository-level pipeline and phase definitions are in the root `Mission.md` and `Architecture.md`.

---

## 1. Objective

- **Scientific:** Reproduce all figures from the DLBCL IMC manuscript "Integrative spatial analysis reveals a hierarchy of cellular organization in diffuse large B-cell lymphoma" (v7.4). 328 treatment-naive tumors, 52 markers across 2 IMC panels (immune + stromal), 12 major cell types -> 30 subpopulations, 5 LME (lymphoma microenvironment) classes.
- **Operational:** Hybrid approach — reuse existing Seurat-converted h5ad objects (47 files on cayuga) for cell typing/labels, consolidate into sc_tools checkpoint format. Reproduce all 13 figures (5 main + 8 supplementary) in a reproducible Snakemake pipeline.
- **Runtime:** cayuga (`/home/fs01/juk4007/elementolab/sc_tools/projects/imc/lymph_dlbcl`). Data already present there.
- **Panels:** T2 (immune) + S2 (stromal) — most complete iterations with DLC_code and cell type labels.

---

## 2. Context and Conventions

- **Remote data (read-only reference):** `/home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2` (cayuga).
- **Local project root:** `./projects/imc/lymph_dlbcl/`.
- **Existing assets:** 47 Seurat-converted h5ad in `results/seurat_converted/`; notebooks in `notebooks/dlbcl_notebooks/DLBCLv2/`; manuscript in `manuscript/`.
- **Target:** Standard checkpoint nomenclature per root `Architecture.md`; **keep immune and stromal panels separate**.
- **Manuscript figures:** `manuscript/Figures_v7.1/` (Adobe Illustrator originals for visual comparison).

---

## 3. Figures to Reproduce

| Figure | Content | Key Data |
|--------|---------|----------|
| **Fig 1** | Single-cell atlas: UMAP, marker heatmap, dotplot, cell proportions | Both panels, celltyped P4 |
| **Fig 2** | 5 LME classes: complex heatmap, abundance, IMC images | TME cluster assignments, cell proportions |
| **Fig 3** | Clinical: KM curves (OS/PFS), COO, mutations | Clinical CSV, mutation table, survival data |
| **Fig 4** | Spatial: 20 communities, neighborhood enrichment | Spatial community h5ad (k=30 kNN) |
| **Fig 5** | ML: ROC, feature importance, IHC validation | TME z-scores, IHC data, clinical |
| **Supp 1** | Panel composition and QC metrics | Marker lists, intensity distributions |
| **Supp 2** | B cell subcluster analysis | B cell subset, re-clustering, marker heatmap |
| **Supp 3** | T cell/myeloid characterization | T cell + myeloid subsets, markers, proportions |
| **Supp 4** | Vessel analysis | Endothelial subsets, vessel density per LME |
| **Supp 5** | TME clustering sensitivity | Alternative k values, stability analysis |
| **Supp 6** | Mutation landscape | Per-gene mutation rates, co-occurrence, LME enrichment |
| **Supp 7** | RNA-protein comparison (CosMx) | IMC protein vs CosMx RNA correlation |
| **Supp 8** | Extended survival and IHC validation | Multivariate Cox, subgroup KM, IHC-based TME prediction |

---

## 4. Implementation Phases

### Phase 0: Data Audit and Organization

- [x] Step 1: Remote file inventory and h5ad conversion (47 files)
- [ ] **0.1** Download missing clinical/metadata from cayuga (`scripts/download_clinical_metadata.sh`)
- [ ] **0.2** Validate key h5ad objects (`scripts/validate_h5ad_objects.py`)
- [ ] **0.3** Map cell type labels to manuscript nomenclature (`metadata/celltype_map_immune.json`, `metadata/celltype_map_stromal.json`)

### Phase 1: AnnData Construction (Per Panel)

- [ ] **1.1** Build panel checkpoints from T2/S2 preprocessing objects (`scripts/build_panel_adata.py`)
- [ ] **1.2** Attach clinical metadata (`scripts/attach_clinical_metadata.py`)
- [ ] **1.3** Build celltyped checkpoints (P4)
- [ ] **1.4** Build spatial community object (`scripts/build_spatial_adata.py`)

### Phase 2: QC and Normalization Verification

- [ ] **2.1** Verify normalization (`scripts/verify_normalization.py`)
- [ ] **2.2** QC report (`scripts/run_qc_report.py`)

### Phase 3: Cell Typing Validation and LME Construction

- [ ] **3.1** Validate cell types (`scripts/validate_celltypes.py`)
- [ ] **3.2** Assign LME classes (`scripts/build_lme_classes.py`)

### Phase 4: Figure Reproduction

- [ ] **4.1** Fig 1: Single-cell atlas (`scripts/fig1_single_cell_atlas.py`)
- [ ] **4.2** Fig 2: LME classes (`scripts/fig2_lme_classes.py`)
- [ ] **4.3** Fig 3: Clinical (`scripts/fig3_clinical.py`)
- [ ] **4.4** Fig 4: Spatial (`scripts/fig4_spatial.py`)
- [ ] **4.5** Fig 5: ML framework (`scripts/fig5_ml_framework.py`)
- [ ] **4.6** Supp 1: QC panels (`scripts/supp_fig1_qc_panels.py`)
- [ ] **4.7** Supp 2: B cell (`scripts/supp_fig2_bcell.py`)
- [ ] **4.8** Supp 3: T cell/myeloid (`scripts/supp_fig3_tcell_myeloid.py`)
- [ ] **4.9** Supp 4: Vessel (`scripts/supp_fig4_vessel.py`)
- [ ] **4.10** Supp 5: TME sensitivity (`scripts/supp_fig5_tme_sensitivity.py`)
- [ ] **4.11** Supp 6: Mutations (`scripts/supp_fig6_mutations.py`)
- [ ] **4.12** Supp 7: RNA-protein (`scripts/supp_fig7_rna_protein.py`)
- [ ] **4.13** Supp 8: Extended survival (`scripts/supp_fig8_extended_survival.py`)

### Phase 5: Snakemake Pipeline

- [ ] **5.1** `config.yaml` with project parameters
- [ ] **5.2** `Snakefile` with full DAG

### Phase 6: Verification

- [ ] **6.1** Visual comparison with `manuscript/Figures_v7.1/` originals
- [ ] **6.2** Numerical validation (cell counts, LME sizes, KM p-values, AUC)
- [ ] **6.3** End-to-end Snakemake run

---

## 5. Completed Tasks

- [x] Create project layout: `data/`, `figures/`, `metadata/`, `scripts/`, `results/`, `outputs/`, and planning docs.
- [x] Step 1: Remote file listing obtained; inventory script run; annotated inventory produced (`metadata/remote_file_inventory_annotated.csv`).
- [x] RDS metadata extracted (47 files + notebook traceability); Seurat-to-h5ad conversion done (48 objects in `results/seurat_converted/`).

---

## 6. Key Files

| What | Where |
|------|-------|
| Immune full (T2) | `results/seurat_converted/tcell_2_preprocessing/t2_SO_seurat.h5ad` (1.63M cells, 49 vars, has DLC_code + labels) |
| Stromal full (S2) | `results/seurat_converted/stroma_2_preprocessing/S1_seurat_SO.h5ad` (1.55M cells, 50 vars) |
| Immune merged | `results/seurat_converted/tcell_merged/SO2_seurat_bcell.h5ad` (979K cells, has DLC_code + labels) |
| Stromal merged | `results/seurat_converted/stroma_merged/SO2_seurat_bcell.h5ad` (1.17M cells, has meta + recluster) |
| Spatial communities | `results/seurat_converted/stroma_spatial/SO_k30_community_cluster.h5ad` (1.65M cells, cluster + community) |
| TME object | `results/seurat_converted/s_tme_seurat.h5ad` (576 obs, 34 vars) |
| h5ad inventory | `metadata/h5ad_inventory_summary.csv` (48 rows) |
| Remote file inventory | `metadata/remote_file_inventory_annotated.csv` |
| Manuscript | `manuscript/DLBCL_hyperion_draft_v7.4_comments.docx` |
| Original figures | `manuscript/Figures_v7.1/` (Adobe Illustrator) |

---

## 7. Risks and Mitigations

1. **Missing sample mapping for stromal panel** — S2 objects lack DLC_code. Mitigation: use cell ID CSVs or orig.ident mapping.
2. **Spatial coordinates** — may be stored differently in Seurat vs AnnData. Mitigation: extract from quant CSVs.
3. **ML reproducibility** — random seed dependent. Mitigation: fix seed, report CI.
4. **R-to-Python porting** — original analysis in R. Mitigation: validate outputs numerically, not pixel-perfect.

---

## 8. Traceability and Reproducibility

- **Scripts:** All automated steps live under `scripts/` with clear inputs and outputs.
- **Manifests:** `metadata/download_manifest_phase0.csv` for downloaded files.
- **Journal:** Log material decisions in `Journal.md`.
- **Snakemake:** Full pipeline in `Snakefile` + `config.yaml`.
