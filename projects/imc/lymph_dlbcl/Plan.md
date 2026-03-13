# Plan: Two-Panel IMC DLBCL (lymph_dlbcl)

## Context

**Scientific objective:** Reproduce all figures from the DLBCL IMC manuscript (v7.4) — 328 treatment-naive tumors, 52 markers across immune (T2) + stromal (S2) panels, 12 major cell types, 30 subpopulations, 5 LME classes.

**Hybrid approach:** Plan A reuses 47 Seurat-converted h5ad objects for cell typing/labels and figure reproduction. Plan B reprocesses all raw Hyperion data from scratch via steinbock on cayuga to produce fully reproducible per-patient quantifications. Both plans run in parallel.

---

## Phase Status

- [x] `ingest_raw` / `ingest_load` — 47 h5ad files converted; panel checkpoints built
- [x] `metadata_attach` — clinical metadata attached (78/138 immune samples matched)
- [x] `celltype_manual` — immune panel celltyped with LME_class (5 classes via Hungarian matching)
- [x] `scoring` — LME classes computed (Cold 31%, Stromal 19%, Cytotoxic 21%, T cell Reg 15%, CD206 13%)
- [x] `biology` — 64 PDFs generated; visual QA pass on all non-blocked figures
- [ ] `qc_filter` — normalization verification and QC report pending
- [ ] Plan B (`ingest_raw` from scratch) — infrastructure ready, steinbock SIF pull in progress

---

## Plan A: Figure Reproduction

| Figure | Status | Notes |
|--------|--------|-------|
| Fig 1: Single-cell atlas | [x] PASS | UMAP x2, heatmap (z-scored 40 markers x 16 subtypes), prevalence; layers['raw'] fallback for X=zeros |
| Fig 2: LME classes | [x] PASS | Heatmap, stacked bar x2, violin, proportions; all 5 classes present |
| Fig 3: Clinical | [x] PASS | KM OS p=0.0147, KM PFS p=0.0024, Cox forest, COO distribution |
| Fig 4: Spatial | [x] PASS | 30 cell types x 18 communities, diversity bar, enrichment bubble, ML ROC+FI |
| Fig 5: CosMx | [x] DONE | 5 panels: UMAP (241K non-B cells), dotplot (24 markers), spatial map, 6 pathways, receptor-ligand |
| Extended: ML framework | [x] PASS | ROC, feature importance, confusion, IHC, metrics (in `extended/`) |
| Supp 1: QC panels | [x] PASS | Markers, distributions, per-sample QC, both panels |
| Supp 2: B cell | [x] PASS | UMAP x4, heatmap x4 |
| Supp 3: T cell/myeloid | [x] PASS | UMAP, heatmap, proportions (T2 + S1 + v2) |
| Supp 4: Vessel | [x] PASS | CD31 violin, vessel density by LME |
| Supp 5: TME sensitivity | [x] PASS | Cluster metrics, stability, k=5/8/10 heatmaps |
| Supp 6: Mutations | [ ] BLOCKED | Needs `CTMA121_mut_table.csv` |
| Supp 7: RNA-protein | [x] DONE | 7 overlapping markers, per-marker scatter + correlation bar (3 common cell types) |
| Supp 8: Extended survival | [x] PASS | Multivariate Cox, KM by COO subgroup |

### Remaining Plan A tasks

- [x] Fig 5 CosMx: 5 panels generated (UMAP, dotplot, spatial, pathways, receptor-ligand)
- [x] Supp 7 RNA-protein: 2 panels generated (per-marker scatter, correlation summary)
- [ ] Obtain `CTMA121_mut_table.csv` for Supp 6 (mutation landscape)
- [ ] Visual comparison with `manuscript/Figures_v7.1/` originals — final pass
- [ ] Numerical validation (cell counts, LME sizes, KM p-values, AUC)
- [ ] End-to-end Snakemake dry-run and full run

---

## Plan B: Raw IMC Reprocessing (steinbock on cayuga)

### Raw data structure

```
/athena/elementolab/scratch/dym2001/data/hyperion/DLBCL/
├── DLBCL_TMA_slide1/          # STROMAL; 11 ROIs; no DLC in filename
├── DLBCL_TMA_slide2/          # IMMUNE; 2 ROIs
├── DLBCL_TMA_slide3/          # STROMAL; duplicates (DLC in filename)
├── DLBCL_TMA_slide4/          # IMMUNE; DLC in filename
├── DLBCL_TMA_slideC_stroma/   # STROMAL; 170 txt; DLC in filename
├── DLBCL_TMA_slideC_tcell/    # IMMUNE; 149 txt; DLC in filename
├── DLBCL_TMA_cornell_stroma/  # STROMAL; 68 txt; well# in filename
├── DLBCL_TMA_cornell_tcell/   # IMMUNE; 81 txt; well# in filename
└── meta/                      # Clinical CSVs, punch notes, mutations
```

**9 steinbock runs** (slide-level, not per-patient). ROI-to-DLC mapping done in post-processing.

### Completed steps

- [x] Panel CSVs: `metadata/panel_immune_t2.csv`, `metadata/panel_stromal_s2.csv` (Metal_Tag + raw_channel_name)
- [x] `config.yaml` `phase0a` block (enabled: false; steinbock_sif/steinbock_img paths)
- [x] `Snakefile` `run_imc_pipeline` rule (steinbock via Apptainer; handles MCD + txt)
- [x] `scripts/setup_steinbock.sh` (pulls SIF on cayuga)
- [x] `scripts/convert_panel_for_steinbock.py` (sc_tools panel CSV to steinbock format)
- [x] Slide-level batch manifests with real cayuga paths
- [x] Files synced to cayuga

### Next steps (require cayuga VPN)

- [ ] Verify steinbock SIF pull completed:
  ```bash
  ssh cayuga "cat ~/elementolab/scratch/steinbock_setup.log | tail -5"
  ```
- [ ] Dry-run:
  ```bash
  cd /home/fs01/juk4007/elementolab/sc_tools/projects/imc/lymph_dlbcl
  snakemake -n phase0a_immune --config 'phase0a={enabled: true}'
  ```
- [ ] Pilot run — `slideC_tcell` (149 ROIs, DLC in filename, simplest mapping):
  ```bash
  snakemake data/processed/slideC_tcell/.done \
    --config 'phase0a={enabled: true}' \
    --executor slurm --default-resources "slurm_partition=scu-cpu" \
    --cores 8 --mem-mb 32000 -j 1
  ```
- [ ] Full run (all 8 slides): `snakemake phase0a_immune phase0a_stromal ...`
- [ ] Write `scripts/aggregate_by_dlc.py` — map ROI names to DLC codes:
  - Slides with DLC in filename: parse `DLC_XXXX` from ROI name
  - Cornell slides: join with `meta/cornell/CTMA 121 punch file (1mm) WT Notes.xlsx`
  - Slide1: join with `meta/BCCA/` mapping
  - Output: `data/{DLC_code}/adata.ingested.h5ad`
- [ ] Enable `phase0a.enabled: true` in config.yaml and submit SLURM jobs

---

## Verification

- [ ] Visual comparison with `manuscript/Figures_v7.1/` originals (final pass)
- [ ] Numerical validation: cell counts, LME sizes, KM p-values, AUC
- [ ] End-to-end Snakemake run (`snakemake --cores 8 all`)
- [ ] Validate h5ad objects (`scripts/validate_h5ad_objects.py`)
- [ ] Map cell type labels to manuscript nomenclature (celltype_map JSONs)
- [ ] Verify normalization (`scripts/verify_normalization.py`)

---

## Blocked Items

- **Supp 6 (mutations):** Needs `CTMA121_mut_table.csv` — check `meta/cornell/` on cayuga
- **Fig 5 + Supp 7 (CosMx):** DONE (2026-03-13). RNA h5ad converted, clusters annotated, 5+2 panels generated
- **Plan B execution:** Steinbock SIF pull in progress on cayuga; all downstream steps blocked until complete
- **`scripts/aggregate_by_dlc.py`:** Not yet written; needed after steinbock runs complete
- **Slide1 DLC mapping:** No DLC code in slide1 filenames; mapping file needed from `meta/BCCA/`

---

## Key Files

| What | Path |
|------|------|
| Immune full (T2) | `results/seurat_converted/tcell_2_preprocessing/t2_SO_seurat.h5ad` (1.63M cells, 49 vars) |
| Stromal full (S2) | `results/seurat_converted/stroma_2_preprocessing/S1_seurat_SO.h5ad` (1.55M cells, 50 vars) |
| Immune merged | `results/seurat_converted/tcell_merged/SO2_seurat_bcell.h5ad` (979K cells) |
| Stromal merged | `results/seurat_converted/stroma_merged/SO2_seurat_bcell.h5ad` (1.17M cells) |
| Spatial communities | `results/seurat_converted/stroma_spatial/SO_k30_community_cluster.h5ad` (1.65M cells) |
| TME object | `results/seurat_converted/s_tme_seurat.h5ad` (576 obs, 34 vars) |
| Clinical metadata | `metadata/DLC380_clinical.tsv` (348 cases) |
| h5ad inventory | `metadata/h5ad_inventory_summary.csv` (48 rows) |
| Manuscript | `manuscript/DLBCL_hyperion_draft_v7.4_comments.docx` |
| Original figures | `manuscript/Figures_v7.1/` (Adobe Illustrator) |
| Figure config | `scripts/figure_config.py` (LME_COLORS, CELLTYPE_COLORS, apply_figure_style) |
| Panel CSVs | `metadata/panel_immune_t2.csv`, `metadata/panel_stromal_s2.csv` |
| Raw data (cayuga) | `/athena/elementolab/scratch/dym2001/data/hyperion/DLBCL/` |
| Project root (cayuga) | `/home/fs01/juk4007/elementolab/sc_tools/projects/imc/lymph_dlbcl` |

---

## Technical Decisions

- **steinbock over imctools+CellProfiler:** steinbock is ElementoLab standard (v0.16.x), handles MCD+txt natively, runs as Apptainer container, produces AnnData directly. Neither imctools nor CellProfiler is installed on cayuga.
- **Cellpose for segmentation:** `latest-cellpose` image tag. Nuclear channels: HistoneH3 (In113Di) + DNA1 (Ir191Di). DeepCell not used (TF/numpy conflict).
- **Slide-level steinbock runs:** 9 runs (not 300+ per-patient) for efficiency; ROI-to-DLC mapping in post-processing.
- **X=zeros in celltyped checkpoint:** Expression data in `layers['raw']` (Seurat scale.data), not X. All 8+ scripts use `layers.get('raw', X)` fallback.
- **Stromal panel patient tracking:** S2 h5ad lacks DLC_code; barcode regex parsing used as fallback.
- **LME assignment:** Hungarian optimal matching on TME h5ad; validates proportions within 5% of manuscript.
- **R-to-Python porting:** Validate numerically, not pixel-perfect. Fix seed=42 for ML reproducibility.
- **Cornell well mapping:** `meta/cornell/CTMA 121 punch file (1mm) WT Notes.xlsx` maps well positions to DLC IDs.

---

## Conventions

- **Panels:** T2 (immune) + S2 (stromal) kept separate; merge only for cross-panel analyses
- **Cell types:** 12 major -> 30 subpopulations; labels from DLC_code (pre-existing)
- **LME classes:** 5 classes per manuscript; stored in `obs['LME_class']`
- **Reference data (read-only):** `/home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2`
- **Figure format:** Match manuscript figure layout exactly for visual QA comparison
- **Statistical correction:** Benjamini-Hochberg (FDR) for all multiple comparisons
