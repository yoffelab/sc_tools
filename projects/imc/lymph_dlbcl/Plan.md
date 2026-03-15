# Plan: Two-Panel IMC DLBCL (lymph_dlbcl)

## Goal

**Scientific objective:** Reproduce all 13 figures (5 main + 8 supplementary) from the DLBCL IMC manuscript "Integrative spatial analysis reveals a hierarchy of cellular organization in diffuse large B-cell lymphoma" (v7.4). 328 treatment-naive tumors, 52 markers across immune (T2) + stromal (S2) panels, 12 major cell types, 30 subpopulations, 5 LME classes.

**Why this matters:** A figure reproduction is only meaningful if the reproduced results confirm (or explicitly challenge) the original biological findings. Code that runs and produces a PDF is necessary but not sufficient. Every figure must demonstrate that the underlying biological signal is present and interpretable.

### 3-Part Figure QA Standard

Every figure must satisfy ALL THREE criteria before it can be marked PASS:

| # | Criterion | What it means | PASS requires |
|---|-----------|---------------|---------------|
| 1 | **Technical** | Code compiles, PDF exists | Script runs without error; expected output files are present |
| 2 | **Biological signal** | Figure shows a clear, interpretable result | A specific, documentable finding (statistically significant difference, meaningful pattern, or an informative null result). Not just "the heatmap rendered" but "the heatmap shows CD8 T cells enriched in Cytotoxic LME with p < 0.01" |
| 3 | **Manuscript alignment** | Reproduced finding is compared against the paper | Reviewer documents whether the result supports, contradicts, or partially matches the manuscript claim. Effect sizes, p-values, and direction-of-effect are explicitly compared |

A figure that is technically correct but has no documented biological signal or no manuscript comparison is **not** PASS. It is "Technical PASS, needs bio review."

---

## Data Construction Status (prerequisite for all figures)

Building properly annotated checkpoint files for both panels is the foundation. Without correct data, figures cannot pass the biological signal criterion.

### Immune Panel (T2) — Functional, partially annotated

| Checkpoint | Exists | Shape | Key Annotations |
|------------|--------|-------|-----------------|
| `adata.immune.raw.p1.h5ad` | YES (3.2 GB) | 1,628,885 x 49 | 138 samples via DLC_code |
| `adata.immune.annotated.p2.h5ad` | YES (3.5 GB) | 1,628,885 x 49 | Clinical on 78/138 samples (39% of cells) |
| `adata.immune.celltyped.p4.h5ad` | YES (3.5 GB) | 1,628,885 x 49 | 19 cell types; LME_class on 40% of cells (287/338 patients) |

**Known issues:**
- X = Seurat ScaleData (centered/scaled, negative values). No raw counts preserved.
- layers['raw'] = copy of X (NOT actual raw counts)
- 60/138 samples lack clinical data (non-DLC IDs like CTMA121_*, I16-*)
- LME_class empty for 60% of cells (samples not in 287-patient assignment table)
- No spatial coordinates (obsm['spatial'] missing; only X_pca from Seurat)
- No UMAP embedding (obsm['X_umap'] missing)

### Stromal Panel (S2) — BROKEN, not usable for per-patient analysis

| Checkpoint | Exists | Shape | Key Annotations |
|------------|--------|-------|-----------------|
| `adata.stromal.raw.p1.h5ad` | YES (3.0 GB) | 1,552,303 x 50 | **sample=0 for ALL cells** (orig.ident not mapped) |
| `adata.stromal.annotated.p2.h5ad` | YES (3.0 GB) | 1,552,303 x 50 | **No clinical data attached** |
| `adata.stromal.celltyped.p4.h5ad` | **MISSING** | - | Never built |
| `adata.stromal.spatial.communities.h5ad` | **MISSING** | - | Never built (source exists in seurat_converted/) |

**Root cause:** `S2_seurat_cellid.csv` was never downloaded from cayuga. This file maps cell indices to DLC patient codes. Without it, all 1.55M stromal cells appear as a single sample. The barcode parsing fallback in `build_panel_adata.py` found no DLC codes in the integer cell indices.

**Cell type status:** All celltype = Unknown. The spatial community object (`SO_k30_community_cluster.h5ad`, 1.65M cells) has 30 subtypes and shares the same integer index format — label transfer is possible but not yet done.

### Steps to fix (all require cayuga VPN)

**Confirmed 2026-03-13:** Exhaustive local search found NO alternative source for cell-to-patient mapping. All stromal h5ad files (full, merged, spatial, subsets) lack DLC_code. The RDS source (`S1_seurat_SO.rds`) also has no DLC_code — it was maintained externally in `S2_seurat_cellid.csv`. There is no overlap-based reconstruction possible.

1. [ ] **Download `S2_seurat_cellid.csv`** from cayuga (`stroma_2_preprocessing/S2_seurat_cellid.csv`, ~60MB) — unlocks stromal patient mapping
2. [ ] Re-run `build_panel_adata.py --panel stromal` with cell ID CSV
3. [ ] Transfer cell type labels from spatial community object (matching index format)
4. [ ] Re-run `attach_clinical_metadata.py --panel stromal`
5. [ ] Build `adata.stromal.celltyped.p4.h5ad` with LME class
6. [ ] Extract spatial coordinates from raw quant CSVs (both panels)
7. [ ] Re-convert RDS with `@assays$Protein@counts` to get raw counts (both panels)

### Impact on figures

Most figures currently use only the immune panel (which works). Figures that need the stromal panel (parts of Fig 1 atlas, Fig 2 composition, spatial analyses) may show incomplete results. The stromal fix is required for a complete manuscript reproduction.

---

## Phase Status

- [x] `ingest_raw` / `ingest_load` -- 47 h5ad files converted; immune panel checkpoints built; **stromal broken (see above)**
- [x] `metadata_attach` -- immune: clinical on 78/138 samples; **stromal: no clinical data**
- [x] `celltype_manual` -- immune: 19 cell types + LME_class (40% coverage); **stromal: all Unknown**
- [x] `scoring` -- LME classes computed (Cold 31%, Stromal 19%, Cytotoxic 21%, T cell Reg 15%, CD206 13%)
- [x] `biology` -- ~66 PDFs generated; first visual QA pass complete (2026-03-09); bio review pending
- [ ] `qc_filter` -- normalization verification and QC report pending
- [ ] Plan B (`ingest_raw` from scratch) -- infrastructure ready, steinbock SIF pull in progress

---

## Plan A: Figure Reproduction

### Figure QA Status

#### Main Figures

| Figure | Technical | Biological Signal | Manuscript Alignment | Status |
|--------|-----------|-------------------|----------------------|--------|
| **Fig 1: Single-cell atlas** | PARTIAL (1b heatmap blank, 1d missing) | Immune UMAP: 20 subtypes (10 T, 6 M, bcell, stroma, other) with distinct clusters. Prevalence: M0_PDL1+ 12.8%, T4_CD4+ 10.2%, M2_CD68+ 9.3%. COO-LME: chi2 p=6.06e-4. **BUT:** 1b heatmap body is blank (z-scoring on already-scaled data → near-zero variance); stromal UMAP shows only 3 categories (Unknown/bcell/stroma); Cold missing from COO plot | Partially supports: 20 immune subtypes plausible but cannot verify full 30 (stromal broken). COO-LME association significant. | **REVISE** — fix heatmap contrast, add fig1d |
| **Fig 2: LME classes** | YES | All 5 LME classes with correct enrichment signatures. Proportions: Cold 31.1% (n=89), Stromal 19.2% (n=55), Cytotoxic 21.3% (n=61), T cell Reg 15.4% (n=44), CD206 12.9% (n=37). n=286 total. Violin: M3_IDO+ 10x enriched in Cytotoxic (p<0.0001); S3_PDPN+ 5x enriched in Stromal (p<0.0001). Heatmap z-scores show Cold=B cell enriched, Cytotoxic=CD8+GrB+/macrophage enriched, CD206=CD206+ macrophage enriched | **Supports:** Rank order preserved. Proportions within ~5pp (largest delta: CD206 12.9% vs ms 8.2%). Enrichment directions all match. n=286 vs ms 328 (42 samples lost in join) | **Near PASS** |
| **Fig 3: Clinical** | YES | OS: log-rank p=0.0147 (n=272). PFS: p=0.0024 (n=272). Cox (univariate, ref=CD206): Cytotoxic HR=1.76 (CI 1.02-3.03, p=0.042); Cold HR=1.20 (ns); Stromal HR=0.79 (ns). COO-LME: chi2 p=1.55e-3. Cytotoxic=most ABC (45%), Stromal=most GCB (75%) | **Partially supports:** PFS closely matches (0.0024 vs ms 0.0025). OS weaker (0.0147 vs ms 0.0064). Cox HR attenuated (1.76 vs ms 2.69) — our model is univariate, ms uses multivariate with age/IPI/stage/COO. Direction correct. | **REVISE** — add covariates to Cox |
| **Fig 4: Spatial** | PARTIAL (4a degraded) | **4a:** Only 3 cell types (Unknown/bcell/stroma) — uninformative due to broken stromal annotations. **4b:** Shannon entropy 0.30-0.79 but on 3 categories. **4c:** Community-LME enrichment shows real signal: Cold=core_tumor enriched, Stromal=C03/C07-C10, Cytotoxic=C06-C09. **4d:** GBM AUC: Cold 0.94, Cytotoxic 0.94, Stromal 0.93, T cell Reg 0.86, CD206 0.80. Top feature: core_tumor fraction (0.155) | **Partially supports:** AUC range 0.80-0.94 matches ms exactly. Community-LME patterns biologically coherent. BUT 4a/4b are uninformative (need stromal fix) | **REVISE** — fix stromal annotations for 4a/4b |
| **Fig 5: CosMx** | YES | 241,471 non-B cells, 5 cell types on UMAP (CD8 T, DC, Endothelial, Fibroblast, Neutrophil). Dotplot: 23 markers x 7 types. Spatial: 78K cells from CTMA121. Pathways: Fibroblasts highest VEGF/TGFb; "Neutrophils" highest IFNg. Receptor-ligand: PD1-PDL1 strongest with Plasma source (0.112); ICOS-ICOSL strongest with B cell source (0.105) | **Partially supports:** Cell count matches. BUT "Neutrophil" cluster expresses CD68/CD163 (likely macrophage); CD4 T and NK missing; UMAP has poor separation (ScaleData in X) | **REVISE** — fix cell type annotations, clean stale ML files |

#### Extended Data

| Figure | Technical | Biological Signal | Manuscript Alignment | Status |
|--------|-----------|-------------------|----------------------|--------|
| **Extended: ML framework** | YES | RF macro AUC=0.90. Per-class: Stromal 0.95, Cytotoxic 0.91, T cell Reg 0.90, CD206 0.87, Cold 0.84. Top features: PDPN (0.115), IDO (0.078), CD163 (0.065), CD4 (0.060). Accuracy=0.664, best recall Stromal (0.84), worst CD206 (0.38). IHC validation is placeholder only | **Partially supports:** AUC 0.84-0.95 within ms range 0.80-0.94. Cold lowest AUC in RF vs highest in GBM (community) — makes biological sense. Feature importance biologically coherent (PDPN=stromal, IDO=immune evasion, CD163=M2) | **Near PASS** (IHC placeholder) |

#### Supplementary Figures

| Figure | Technical | Biological Signal | Manuscript Alignment | Status |
|--------|-----------|-------------------|----------------------|--------|
| **Supp 1: QC panels** | YES (6/6 PDFs) | Panel (a): mean marker intensity bars for 49-50 markers per panel. Panel (b): per-marker intensity histograms showing bimodal vs continuous distributions. Panel (c): per-sample nCount/nFeature boxplots. Script uses `layers['raw']` fallback | Standard IMC QC — matches manuscript expectations. Stromal per-sample QC (1c) uninformative due to sample-0 bug (single box). ScaleData risk if `layers['raw']` absent | **Near PASS** (stromal QC degraded) |
| **Supp 2: B cell** | YES (8/8 PDFs) | UMAP + z-scored marker heatmap for B cell subclusters across 4 objects (immune T2, immune merged, stromal S2, stromal merged). Reveals molecular heterogeneity within tumor B cells (proliferation, activation, GCB/ABC markers) | Matches manuscript — B cell subclustering central to DLBCL. Risk of double z-scoring (z-scoring already-scaled ScaleData). Stromal objects inherit sample-0 bug | **Near PASS** (double z-score risk) |
| **Supp 3: T cell/myeloid** | YES (12/12 PDFs) | UMAP + heatmap + per-sample proportions for 4 subsets (tcell T2, tcell T2 v2, myeloid T2, myeloid S1). T cell subsets separate CD4/CD8/Treg/exhausted/effector. Myeloid separates M1/M2/DC. Per-sample proportions link to LME classes | Good alignment — T cell and myeloid subtyping central to TME hierarchy. myeloid_S1 proportions broken (stromal sample bug). Same double z-scoring risk | **Near PASS** (myeloid_S1 degraded) |
| **Supp 4: Vessel** | PARTIAL (1/2 PDFs) | Panel (a): CD31 violin across cell types — validates endothelial annotation. **Panel (b) MISSING:** vessel density per LME class not generated (LME mapping likely failed — check `lme_class_assignments.csv` join) | Partially supports — endothelial markers validated. Key analysis (vessel density by LME) missing. Script combines annotation-based + CD31-high detection | **REVISE** — debug supp4b generation |
| **Supp 5: TME sensitivity** | YES | Silhouette scores weak (0.06-0.13) with no clear peak at k=10. Bootstrap stability moderate (ARI ~0.50-0.53 across k=5/8/10). k=5 heatmap maps to 5 LME classes; k=10 shows biologically coherent splits (B cell, stromal, myeloid, T cell subtypes) | **Partially supports:** Provides transparency about clustering uncertainty. k=10 not strongly justified by either metric alone — manuscript should note domain knowledge drove the choice. k=5 heatmap recovers the 5 LME classes | **REVISE** — acknowledge weak silhouette |
| **Supp 6: Mutations** | NO -- BLOCKED | N/A | N/A | BLOCKED (needs CTMA121_mut_table.csv) |
| **Supp 7: RNA-protein** | YES | 7 overlapping markers, 3 matched cell types (B cell, CD8 T, NESC). 6/7 markers show r>0.91 (CD68 r=1.00 p=0.036; MS4A1 r=1.00 p=0.053; CD4 r=0.99; PECAM1 r=0.95; ICOS r=0.91). **FOXP3 r=-1.00 (p=0.007)** — anti-correlation likely caused by Seurat ScaleData in X (negative centered values flip correlation sign) | **Partially supports:** Positive correlations for lineage markers expected and confirmed. FOXP3 anti-correlation is a red flag — needs fix (use raw counts or acknowledge ScaleData artifact). n=3 cell types critically underpowered (expected 8) | **REVISE** — fix ScaleData issue, improve cell type matching |
| **Supp 8: Extended survival** | YES | Multivariate Cox (ref=CD206): Cytotoxic HR=1.95 (CI 1.07-3.54, p=0.028); Cold HR=1.40 (ns); Stromal HR=1.07 (ns); T cell Reg HR=1.16 (ns). Age HR=1.04 (p=1.3e-5), IPI HR=1.33 (p=1.5e-4), LDH HR=1.0004 (p=0.005). KM by COO: ABC n=87 (no clear separation), GCB n=153 (Cytotoxic worst), U n=26 (underpowered) | **Partially supports:** Cytotoxic worst prognosis confirmed (p=0.028). Clinical covariates independently prognostic. COO subgroups underpowered. Note: IPI already includes age/LDH → possible collinearity. Cytotoxic=worst OS is counterintuitive but consistent with manuscript | **REVISE** — add logrank p to KM plots, fix collinearity |

### Summary (updated 2026-03-13)

- **Technical PASS:** 11/13 (Fig 1 partial — heatmap blank/1d missing; Supp 6 blocked)
- **Biological signal documented:** 12/12 non-blocked figures now have specific quantitative findings
- **Manuscript alignment verified:** 12/12 non-blocked figures compared (2 Near PASS, 10 Partially supports)
- **Near PASS (all 3 criteria):** 2/13 (Fig 2, Extended ML)
- **REVISE:** 9/13 (specific fixes identified per figure)
- **BLOCKED:** 1/13 (Supp 6 — mutation data)
- **Near PASS:** 5/13 (Fig 2, Extended ML, Supp 1, Supp 2, Supp 3)
- **Note:** Supp 1-3 "Near PASS" contingent on ScaleData/double-z-scoring not corrupting results

### Priority Fixes (highest impact)

1. **Fix Fig 1b heatmap** — z-scoring on Seurat ScaleData produces blank heatmap. Use `layers['raw']` with robust scaling or plot raw ScaleData without additional z-scoring
2. **Fix CosMx cell type annotations (Fig 5)** — "Neutrophil" cluster expresses CD68/CD163 (macrophage markers); CD4 T and NK missing. Re-run `annotate_rna_clusters.py` with refined marker sets
3. **Add covariates to Cox model (Fig 3c)** — current univariate HR=1.76 vs ms multivariate HR=2.69. Add age, IPI, stage, COO
4. **Fix Supp 7 FOXP3 anti-correlation** — ScaleData in X (negative centered values) flips sign. Need raw counts or explicit handling
5. **Download `S2_seurat_cellid.csv` from cayuga** — unlocks stromal panel (Fig 1a stromal, Fig 4a/4b, full 30 cell types)
6. **Clean stale PDFs** — old ML files in fig5/, duplicate fig2b/fig2d, supp7 placeholder

### Remaining Plan A Tasks

- [x] Fig 5 CosMx: 5 panels generated (2026-03-13)
- [x] Supp 7 RNA-protein: 2 panels generated (2026-03-13)
- [x] **Systematic biological signal review** — completed 2026-03-13 (4 parallel reviewer agents)
- [x] **Manuscript alignment review** — completed 2026-03-13 (documented per-figure comparisons)
- [ ] Obtain `CTMA121_mut_table.csv` for Supp 6 (mutation landscape)
- [ ] Fix Fig 1b heatmap rendering (blank body)
- [ ] Fix CosMx cell type annotations (Neutrophil → Macrophage, add CD4 T/NK)
- [ ] Add multivariate covariates to Fig 3c Cox model
- [ ] Fix Supp 7 ScaleData artifact (FOXP3 anti-correlation)
- [ ] Add logrank p-values to Supp 8 KM plots; fix IPI/age/LDH collinearity
- [ ] Download `S2_seurat_cellid.csv` from cayuga (VPN required)
- [ ] Clean up stale PDFs
- [ ] End-to-end Snakemake dry-run and full run
- [ ] Final visual QA pass

---

## Figure Verification Protocol

### How to evaluate each criterion

**1. Technical (binary YES/NO):**
- Run the script (or `snakemake <target>`). Does it complete without error?
- Do the expected PDF files exist in `figures/manuscript/`?
- If YES to both, mark Technical = YES.

**2. Biological signal (requires reviewer judgment):**
- Open each PDF and describe the key finding in one sentence.
- The finding must be specific and quantitative where possible: p-values, effect sizes, counts, directions.
- Examples of PASS: "CD8 T cells are 3x enriched in Cytotoxic LME (p < 0.001, BH-corrected)"
- Examples of FAIL: "The heatmap rendered correctly" (no biological content), "Colors look right" (subjective, no finding)
- A meaningful null result is acceptable if documented: "No significant difference in vessel density between Cold and Stromal LME (p = 0.42) -- suggests vessel remodeling is not LME-specific"
- Document the finding in the "Biological Signal" column of the table above.

**3. Manuscript alignment (requires comparison to original):**
- For each documented biological finding, look up the corresponding claim in the manuscript (v7.4) or original figure (`manuscript/Figures_v7.1/`).
- Document one of:
  - **Supports:** Reproduced finding matches manuscript claim (same direction, comparable magnitude)
  - **Partially supports:** Same direction but different effect size, or significant vs non-significant
  - **Contradicts:** Opposite direction or fundamentally different pattern
  - **Cannot compare:** Manuscript claim unclear or data not directly comparable
- For quantitative results, report both values: "Our OS p=0.0147 vs manuscript p=0.02 (same direction, comparable)"
- A contradiction is not automatically a FAIL -- it may reveal a genuine finding. But it must be documented and explained.

### What constitutes PASS vs FAIL

| Scenario | Verdict |
|----------|---------|
| All 3 criteria satisfied | **PASS** |
| Technical YES, bio signal documented, alignment shows contradiction with explanation | **PASS** (with note) |
| Technical YES, bio signal present but not yet documented | **Technical PASS, needs bio review** |
| Technical YES, bio signal documented, alignment not yet checked | **Technical + Bio PASS, needs alignment review** |
| Technical YES, but no interpretable biological signal in the figure | **FAIL** (figure needs revision) |
| Technical YES, bio signal contradicts manuscript with no explanation | **FAIL** (needs investigation) |
| Technical NO (script errors or no PDF) | **FAIL** |

### Process

1. Assessments are documented inline in the figure QA table above (not in separate files).
2. Each review pass should be dated in the Journal.md entry.
3. A figure moves from "Technical PASS" to full "PASS" only after a reviewer fills in both the Biological Signal and Manuscript Alignment columns with specific, substantive content.
4. The first visual QA pass (2026-03-09) covered technical correctness. The biological signal and manuscript alignment review is the next required step.

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
- [ ] Pilot run -- `slideC_tcell` (149 ROIs, DLC in filename, simplest mapping):
  ```bash
  snakemake data/processed/slideC_tcell/.done \
    --config 'phase0a={enabled: true}' \
    --executor slurm --default-resources "slurm_partition=scu-cpu" \
    --cores 8 --mem-mb 32000 -j 1
  ```
- [ ] Full run (all 8 slides): `snakemake phase0a_immune phase0a_stromal ...`
- [ ] Write `scripts/aggregate_by_dlc.py` -- map ROI names to DLC codes:
  - Slides with DLC in filename: parse `DLC_XXXX` from ROI name
  - Cornell slides: join with `meta/cornell/CTMA 121 punch file (1mm) WT Notes.xlsx`
  - Slide1: join with `meta/BCCA/` mapping
  - Output: `data/{DLC_code}/adata.ingested.h5ad`
- [ ] Enable `phase0a.enabled: true` in config.yaml and submit SLURM jobs

---

## Blocked Items

- **Supp 6 (mutations):** Needs `CTMA121_mut_table.csv` -- check `meta/cornell/` on cayuga
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
| CosMx RNA | `data/cosmx/cosmx_rna.h5ad` (464K cells, 1000 genes) |
| CosMx Protein | `data/cosmx/cosmx_prt.h5ad` (687K cells, 67 markers) |
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
