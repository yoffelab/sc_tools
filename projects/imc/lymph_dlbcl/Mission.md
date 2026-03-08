# Mission: Two-Panel IMC DLBCL (lymph_dlbcl)

**Project:** `projects/imc/lymph_dlbcl`
**Current Status:** Manuscript reproduction — 64 PDFs generated, visual QA in progress
**Author:** Junbum Kim
**Last Updated:** 2026-03-06

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
| **Fig 5** | CosMx spatial transcriptomics (deferred) | CosMx data TBD |
| **Ext Data** | ML: ROC, feature importance, IHC validation | TME z-scores, IHC data, clinical |
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
- [x] **0.1** Clinical metadata: `metadata/DLC380_clinical.tsv` (348 cases, direct use — no download dependency)
- [ ] **0.2** Validate key h5ad objects (`scripts/validate_h5ad_objects.py`)
- [ ] **0.3** Map cell type labels to manuscript nomenclature (`metadata/celltype_map_immune.json`, `metadata/celltype_map_stromal.json`)

### Phase 1: AnnData Construction (Per Panel)

- [x] **1.1** Build panel checkpoints from T2/S2 preprocessing objects (`scripts/build_panel_adata.py`)
  - Immune: 1,628,885 cells x 49 markers, 138 samples, 19 cell types
  - Stromal: 1,552,303 cells x 50 markers — **BLOCKER: only 1 sample** (orig.ident=0; needs S2_seurat_cellid.csv)
- [x] **1.2** Attach clinical metadata (`scripts/attach_clinical_metadata.py`) — 78/138 immune samples matched (638K cells)
- [x] **1.3** Build celltyped checkpoints (P4) — immune panel complete with LME_class
- [x] **1.4** Build spatial community object (`scripts/build_spatial_adata.py`)

### Phase 2: QC and Normalization Verification

- [ ] **2.1** Verify normalization (`scripts/verify_normalization.py`)
- [ ] **2.2** QC report (`scripts/run_qc_report.py`)

### Phase 3: Cell Typing Validation and LME Construction

- [ ] **3.1** Validate cell types (`scripts/validate_celltypes.py`)
- [x] **3.2** Assign LME classes (`scripts/build_lme_classes.py`) — 5 classes via Hungarian optimal matching on TME h5ad
  - Cold 31.0% (ms: 35.1%), Stromal 19.2% (21.3%), Cytotoxic 21.3% (20.7%), T cell Reg 15.3% (14.6%), CD206 13.2% (8.2%)

### Phase 4: Figure Reproduction (64 PDFs generated, visual QA in progress)

- [x] **4.1** Fig 1: Single-cell atlas — 4 panels (UMAP x2, heatmap, prevalence; fig1d B cell markers skipped). **v9 fix**: layers['raw'] fallback for X=zeros in p4; heatmap now shows z-scored 40 markers x 16 T/M subtypes with clear biological patterns
- [x] **4.2** Fig 2: LME classes — 6 panels (heatmap, stacked bar x2, violin, proportions x2)
- [x] **4.3** Fig 3: Clinical — 5 panels (KM OS p=0.0147, KM PFS p=0.0024, Cox forest, COO x2)
- [x] **4.4** Fig 4: Spatial — 5 panels (community composition, diversity, enrichment, ML, by LME)
- [x] **4.5** Fig 5: ML framework — 5 panels (ROC, feature importance, confusion, IHC, metrics) → extended/
- [x] **4.6** Supp 1: QC panels — 6 panels (markers, distributions, per-sample QC, both panels)
- [x] **4.7** Supp 2: B cell — 8 panels (UMAP x4, heatmap x4)
- [x] **4.8** Supp 3: T cell/myeloid — 12 panels (UMAP, heatmap, proportions, T2 + S1 + v2)
- [x] **4.9** Supp 4: Vessel — 1 panel (endothelial markers)
- [x] **4.10** Supp 5: TME sensitivity — 5 panels (cluster metrics, stability, k=5/8/10 heatmaps)
- [ ] **4.11** Supp 6: Mutations — **BLOCKED** (needs CTMA121_mut_table.csv)
- [x] **4.12** Supp 7: RNA-protein — placeholder
- [x] **4.13** Supp 8: Extended survival — multivariate Cox, subgroup KM

### Phase 5: Snakemake Pipeline

- [x] **5.1** `config.yaml` with project parameters
- [x] **5.2** `Snakefile` with full DAG

### Plan B: Raw IMC Reprocessing (Phase 0a) — infrastructure ready, blocked on cayuga VPN

- [x] **B.1** `scripts/generate_panel_csv.py` — extracts markers from h5ad var_names → `metadata/panel_immune_t2.csv`, `metadata/panel_stromal_s2.csv` (Metal_Tag TBD from MCD headers)
- [x] **B.2** `scripts/generate_phase0_manifests.py` — builds `metadata/phase0/batch1_immune.tsv` (84 samples), `batch1_stromal.tsv` (349 from clinical fallback); mcd_file paths are PLACEHOLDERS
- [x] **B.3** `config.yaml` `phase0a` block added (enabled: false by default)
- [x] **B.4** `Snakefile` `phase0a_immune`, `phase0a_stromal`, `run_imc_pipeline` rules (gated by phase0a.enabled)
- [ ] **B.5 BLOCKED** Cayuga VPN required to: (a) verify MCD file paths, (b) check imctools/CellProfiler availability, (c) fill Metal_Tag in panel CSVs from MCD headers, (d) run dry-run, (e) submit SLURM jobs

### Phase 6: Verification

- [ ] **6.1** Visual comparison with `manuscript/Figures_v7.1/` originals — QA in progress. Critical fix applied: X=zeros in p4 checkpoint resolved by layers['raw'] fallback in 8 scripts. Fig 1 heatmap/UMAP/prevalence verified. Supp fig 1-2 verified.
- [ ] **6.2** Numerical validation (cell counts, LME sizes, KM p-values, AUC)
- [ ] **6.3** End-to-end Snakemake run

---

## 5. Completed Tasks

- [x] Create project layout: `data/`, `figures/`, `metadata/`, `scripts/`, `results/`, `outputs/`, and planning docs.
- [x] Step 1: Remote file listing obtained; inventory script run; annotated inventory produced (`metadata/remote_file_inventory_annotated.csv`).
- [x] RDS metadata extracted (47 files + notebook traceability); Seurat-to-h5ad conversion done (48 objects in `results/seurat_converted/`).
- [x] **v2 overhaul (2026-03-06):** Figure quality + data pipeline fixes.
  - Added `scripts/figure_config.py` — shared LME_COLORS, LME_ORDER, CELLTYPE_COLORS, COO_COLORS, apply_figure_style().
  - Rewrote `attach_clinical_metadata.py` — uses `metadata/DLC380_clinical.tsv` directly (no download dependency); standardized column names (OS_time, PFS_time, COO, etc.).
  - Rewrote `build_lme_classes.py` — tries patient_tme CSV, then cluster mapping, then k-means fallback; validates all 5 LME classes with manuscript proportions.
  - Rewrote Fig 1-4 scripts with z-scored heatmaps, diverging colormaps, consistent LME/celltype colors, BH-corrected statistics, KM with log-rank p-values, Cox forest with HR/CI.
  - Moved Fig 5 (ML) to `figures/manuscript/extended/ml_classification/` (manuscript Fig 5 is CosMx, deferred).
  - Updated all supp fig scripts to use figure_config.
  - Added `build_panel_adata.py` barcode parsing fallback for stromal panel (DLC_code from cell barcode regex).
  - Updated `config.yaml` to point to `metadata/DLC380_clinical.tsv`.
  - Added skills.md Section 12.8: Figure Intent, Insight, and Readability (general + reproduction-specific).
  - All scripts lint clean (ruff).

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
