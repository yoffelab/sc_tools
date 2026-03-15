---
status: active
created: 2026-03-12
summary: "IMC Lung Cancer Atlas manuscript masterplan — pan-lung spatial proteomics atlas combining ElementoLab + public datasets, 5 main figures + 8-10 supplementary, targeting Nature/Cancer Cell"
---

# IMC Lung Cancer Atlas — Storyboard & Figure Plan

## Context

**Goal:** Create a pan-lung cancer spatial proteomics atlas paper combining our internal ElementoLab IMC data with publicly available datasets. The paper tells two intertwined stories: (1) disease progression from healthy lung through pre-invasive GGO to invasive cancer, and (2) comprehensive TME taxonomy and archetype classification across lung cancer subtypes.

**Working Title:** "A Spatial Proteomics Atlas of the Lung Cancer Microenvironment Reveals Tissue Architecture Remodeling Across Disease Progression and Histological Subtypes"

**Target:** Nature / Cancer Cell | 5 main figures + 8-10 supplementary

**Project directory:** `projects/imc/imc_luad_atlas/`

---

## Data Inventory

### Internal (ElementoLab)

| Dataset | Cancer type | N patients | Markers | Status | Key assets |
|---------|------------|-----------|---------|--------|------------|
| **ggo-imc** | GGO/early LUAD (AAH, AIS, MIA, IAC) | ~40 | Panel G (~35, immune/checkpoint) + Panel H (~35, innate/epithelial) | Cell-typed, UTAG domains, patient risk groups, figures done | Clinical (pathology, radiology, mutations EGFR/KRAS/TP53), EMT analysis, cytokine profiling |
| **ggo_human** | GGO (legacy) | ~10 | ~35 | Cell-typed (legacy naming) | steinbock pipeline, MCD files, segmentation |
| **healthy_lung** | Healthy + COPD reference | 6 (4 healthy, 2 COPD) | 28 (AQP1, aSMA, CD4, CD8A, CD68, CD163, SFTPC, etc.) | UTAG domains, messenger-passed | Proximal/distal location, age groups |

### Public — High Priority (single-cell data available or likely)

| Dataset | Citation | Cancer type | N patients | Markers | Data access | Download | Status |
|---------|----------|------------|-----------|---------|-------------|----------|--------|
| **Sorin et al.** | Nature 2023 | LUAD (5 histological patterns) | 416 | 19 (not 35+) | Public | [Zenodo 7760826](https://zenodo.org/records/7760826) — LungData.zip (2.1 GB) | [x] Downloaded |
| **Desharnais et al.** | Nat Commun 2025 | NSCLC (51 LUAD + 51 LUSC) | 102 | 19 (same Walsh lab panel) | Public | [Zenodo 14625562](https://zenodo.org/records/14625562) — LUSC_LUADzenodo.zip (740 MB) | [x] Downloaded |
| **Hartner et al.** | Nat Commun 2025 | LUAD (driver mutations) | 157 (subset of Sorin 416) | 19 (same panel) | Public | Same as Sorin (Zenodo 7760826) + mutation annotations in supplementary | [x] Included in Sorin download |
| **Cords et al. (CAF)** | Cancer Cell 2024 | NSCLC | 1,070 | 45 | Public | [Zenodo 7961844](https://zenodo.org/records/7961844) v1 — SCE RDS + masks (~13 GB). v3 (15399046) only has mask corrections. | [x] Downloaded (13 GB, SCE format, Zenodo 7961844 v1) |
| **Zhu et al.** | Cancer Cell 2025 | LUAD precursors (TIM-3 focus) | 114 patients, 1618 ROIs | 34 | Synapse (MTA required) | [syn61802885](https://www.synapse.org/Synapse:syn61802885) — requires access request + MTA from MD Anderson | [ ] Blocked — MTA required, skipping |
| **TRACERx IMC** | Cancer Discov 2024 | NSCLC (early-stage) | 81 (198 samples) | 2 panels, 2.3M cells | Restricted (EGA DAC) | Apply: EGAC00001000632; PHLEX test data public on [Zenodo 7973724](https://zenodo.org/records/7973724) | [ ] Not started |
| **CODEX NSCLC** | J Transl Med 2024 | NSCLC (immunotherapy) | ~20+ | 36 (CODEX) | Public | [Zenodo 10258578](https://zenodo.org/records/10258578) — h5ad + 42 ome.tiff (50 GB) | [x] Downloaded |
| **MultiModal NSCLC** | Synapse | NSCLC (ICI response) | TBD | IMC + radiology + genomics | Public (Synapse) | [syn26642505](https://www.synapse.org/Synapse:syn26642505) | [ ] Not started |

**Key discovery:** Sorin, Desharnais, and Hartner are all from the Walsh lab (McGill) with the same 17-marker panel (17 protein markers + DNA channels — originally reported as 19). Hartner's 157 patients are a subset of Sorin's 416 with added mutation data. Total unique from Walsh lab: ~467 patients (416 LUAD + 51 LUSC).

**Cords data format:** R SingleCellExperiment (SCE) in .rds files — needs conversion to AnnData via `anndata2ri` or `zellkonverter`.

### Data Format Details

| Dataset | Raw Format | Key Contents | Cells | ROIs | Conversion Needed |
|---------|-----------|--------------|-------|------|--------------------|
| Sorin 2023 | MATLAB `.mat` + `.tif` masks | `LUAD_IMC_Segmentation/`, `LUAD_IMC_CellType/`, `LUAD_IMC_MaskTif/`, clinical `.xlsx` | 2,141,855 | 536 | DONE |
| Desharnais 2025 | MATLAB `.mat` + `.tif` masks | Same pipeline as Sorin: `segmentation/`, `cellType/`, `LUADTif/`, `LUSCTif/`, `MarkerOrders.xlsx` | 692,746 | 207 (104 LUAD + 103 LUSC) | DONE |
| Cords 2024 | R `.rds` (SingleCellExperiment) | `merge_all_TMAs_sce_RAW.rds` (4.3 GB), `sce_all_annotated.rds` (9.6 GB), metadata CSVs | 5,984,454 (annotated) / 5,984,703 (raw) | Multiple TMAs | DONE |
| COVID Lung IMC | **h5ad** (ready) | `results/covid-imc.h5ad`: 664K cells x 43 features, 36 protein markers | 664,006 (97,098 healthy subset) | 237 / 23 samples | DONE (healthy subset extracted) |
| CODEX NSCLC | **h5ad** (ready) | `anndata_4301_v4_annotated.h5ad`: 211K cells x 27 features, 42 OME-TIFFs | 210,945 | 42 regions | NONE — but CODEX not IMC, partial marker overlap |

### Converted Data Inventory

All converted files at `/athena/elementolab/scratch/juk4007/imc_atlas_data/<dataset>/`

| Dataset | File | Size | Cells | Features | Key Markers | obsm | Key obs columns |
|---------|------|------|-------|----------|-------------|------|-----------------|
| Sorin 2023 | `adata.sorin.h5ad` | 254 MB | 2,141,855 | 17 | TTF1, PanCK, CD68, CD3, CD4, FoxP3, CD8a, CD163, CD20, MPO, CD16, CD14, CD31, CD117, CD94 | spatial | sample_id, cell_type, dataset |
| Desharnais 2025 | `adata.desharnais.h5ad` | 83 MB | 692,746 | 17 | Same 17-marker panel as Sorin | spatial | sample_id, cell_type, cancer_type, dataset |
| Cords 2024 (annotated) | `adata.cords_annotated.h5ad` | 12 GB | 5,984,454 | 43 | MPO, SMA, FAP, HLA-DR, CD146, CD20, CD68, CD3, IDO, Collagen I+Fibronectin, etc. | none (coords in obs Center_X/Center_Y) | ImageNumber, TmaID, Panel, BatchID, CellID |
| Cords 2024 (raw) | `adata.cords_raw.h5ad` | 2.9 GB | 5,984,703 | 43 | Same as annotated | none | Same as annotated |
| COVID healthy | `covid-imc-healthy.h5ad` | — | 97,098 | 43 | 36 protein markers (10-13/19 overlap with Walsh panel) | spatial, X_pca, X_umap | sample, disease, roi, metacluster_label |

**Total cells across converted datasets: ~8.9M**

**Notes on converted data:**
- Sorin+Desharnais share a 17-marker panel (not 19 as originally stated — 17 protein markers + DNA channels)
- Cords spatial coordinates are in obs columns (Center_X, Center_Y), not obsm — will need to move to obsm['spatial'] during harmonization
- COVID healthy marker overlap with Walsh panel — exact matches: CD3, CD4, CD11c, CD20, CD31, CD45, CD68, CD163, Ki67, Vimentin; near-matches: CD8a~CD8, Keratin818~PanCK, AlphaSMA~SMA

### Public — Secondary (SCLC + reference)

| Dataset | Citation | Notes |
|---------|----------|-------|
| **SCLC IMC** | BMC 2025 | 221 ROIs, 37 markers, mouse models — include for pan-lung SCLC comparison |
| **SCLC spatial proteomics** | 2025 | Human SCLC, multiregional (center, margin, peri-tumor), survival associations |
| **Zhu mouse precancers** | Adv Sci 2026 | 284 IMC ROIs, 1.4M cells, mouse LUAD progression models |
| **COVID lung IMC** | Nature 2021 | ElementoLab, 36 markers, healthy lung reference, [Zenodo 4139443](https://zenodo.org/records/4139443) — [x] Downloaded (978 MB, 238 TIFF mask files in 23 sample directories) |
| **ICI spatial signatures** | Nat Genet 2025 | 234 NSCLC, spatial proteomics (n=67) + spatial tx (n=131), response/resistance sigs |
| **LuCA scRNA atlas** | cellxgene | 1.2M cells, 309 patients — reference for cell type annotation transfer |
| **HTAN lung** | Synapse | Multi-modality lung cancer data; check for IMC-specific datasets |

### Data Download Strategy

All downloads at HPC: `/athena/elementolab/scratch/juk4007/imc_atlas_data/`

**Completed downloads:**
- `sorin_2023/LungData.zip` — 2.1 GB (Zenodo 7760826)
- `desharnais_2025/LUSC_LUADzenodo.zip` — 740 MB (Zenodo 14625562)
- `covid_lung_imc/` — 978 MB, 238 TIFFs (Zenodo 4139443)
- `codex_nsclc/` — 50 GB, h5ad + 42 ome.tiff (Zenodo 10258578)
- `cords_2024/` — 13 GB, SCE RDS + metadata CSVs + corrected masks + documentation (Zenodo 7961844 v1)

**Blocked / not started:**
- Zhu et al. (Synapse syn61802885) — requires MTA from MD Anderson, skipping
- TRACERx (EGA) — requires DAC application
- MultiModal NSCLC (Synapse) — not started

### Estimated Total Scale

If we secure key datasets: **~2,000+ patients, ~5M+ cells, pan-lung (LUAD, LUSC, SCLC, precursors, healthy)**

---

## Story Arc

The paper builds from **atlas construction** → **cell phenotypes** → **spatial organization** → **disease progression** → **clinical translation**.

Our unique contributions:
1. **Largest integrated IMC lung cancer atlas** (multi-cohort, multi-center)
2. **Pre-invasive to invasive continuum** (GGO data — unique angle no public dataset covers well)
3. **Pan-lung comparison** (LUAD vs LUSC vs SCLC TME differences)
4. **Healthy lung baseline** (enables progression trajectory)
5. **Spatial archetype taxonomy** with clinical associations

---

## Figure 1: Atlas Construction and Multi-Cohort Integration

**Question:** Can we build a unified spatial proteomics reference across multiple IMC cohorts with partially overlapping panels?

| Panel | Content | Viz type | Data source |
|-------|---------|----------|-------------|
| **a** | Study design schematic: cohorts, patient counts, disease stages covered, total cells | Graphical (BioRender/Illustrator) | All — metadata |
| **b** | Marker panel overlap across all cohorts (shared vs unique markers) | UpSet plot or tile heatmap (markers x datasets) | All panel CSVs |
| **c** | Integrated UMAP — left: colored by cell type, right: colored by cohort | UMAP (2 panels) | All cohorts, shared marker subspace |
| **d** | Marker x cell type expression validation heatmap | Clustered heatmap (z-scored) | Integrated atlas |
| **e** | Cell type proportions by cohort and disease stage | Stacked bar chart | All cohorts |

**Methods:** Harmony on arcsinh-transformed shared markers (~12-15 core: CD3, CD4, CD8, CD68, CD163, PanCK, Ki67, aSMA, CD31, CD20, FoxP3, HLADR). Integration benchmarked with kBET, LISI, silhouette. For non-overlapping markers: UINMF or StabMap for partial-overlap integration/imputation.

**Stats:** Silhouette score (batch vs celltype), chi-square for composition differences across cohorts.

---

## Figure 2: Cell Phenotype Landscape

**Question:** What is the full taxonomy of cell states in the lung cancer TME, and how do they differ across cancer subtypes?

| Panel | Content | Viz type | Data source |
|-------|---------|----------|-------------|
| **a** | Fine-grained cell state UMAP with sub-clustering (T cells, myeloid, epithelial, stromal insets) | Large UMAP + zoom insets | Integrated atlas |
| **b** | T cell exhaustion-activation spectrum: functional markers across T cell subtypes | Dot/bubble plot (subtypes x markers: GrzB, PD1, TIM3, LAG3, CTLA4) | ggo-imc Panel G, Sorin, TRACERx |
| **c** | Myeloid polarization: M1 vs M2 score by disease stage | Density scatter + marginal histograms | All cohorts with CD68/CD163/CD206 |
| **d** | CAF phenotype diversity (if 1070-pt dataset available) or stromal subtype taxonomy | Heatmap of stromal subtypes x functional markers | CAF study (Cancer Cell 2024), ggo-imc |
| **e** | Cross-cohort cell state concordance | Sankey diagram (per-cohort labels → integrated labels) | All cohorts |

**Methods:** Leiden sub-clustering within major lineages. Marker-based annotation. Label transfer via kNN in shared embedding.

**Stats:** ARI, NMI for label concordance. Kruskal-Wallis across disease stages per marker.

---

## Figure 3: Spatial Tissue Architecture and Microenvironment Niches

**Question:** How is the spatial organization of the lung tumor TME structured, and what tissue archetypes emerge?

| Panel | Content | Viz type | Data source |
|-------|---------|----------|-------------|
| **a** | Representative IMC images across disease stages (healthy → GGO → LUAD → LUSC) with cell type overlays | Multi-channel composites (PanCK=R, CD3=G, DNA=B) + cell type maps | ggo-imc, healthy_lung, public |
| **b** | Spatial domain taxonomy (UTAG micro-environments) | Clustered heatmap (domains x cell type fractions) + example ROI domain maps | ggo-imc UTAG, healthy_lung UTAG |
| **c** | Cell-cell neighborhood enrichment | Co-localization z-score heatmap + spatial example plots | All cohorts with spatial coords |
| **d** | Tissue archetype classification (4-6 archetypes: immune-excluded, immune-infiltrated, TLS-rich, myeloid-enriched, stromal-dominated, immune-desert) | Archetype composition heatmap + prevalence stacked bar by disease stage | All cohorts (patient-level) |
| **e** | Spatial diversity metrics (Shannon entropy, Simpson) across stages | Violin plots per disease stage | All cohorts with spatial coords |

**Methods:** UTAG (Bao et al.) for spatial domains. Squidpy for neighborhood enrichment (radius=40um). NMF/consensus clustering for archetype discovery (analogous to LME classes in lymph_dlbcl). Permutation tests for enrichment.

**Stats:** Permutation test with FDR for neighborhood enrichment. Fisher exact for archetype-stage associations. Kruskal-Wallis + Dunn post-hoc for diversity.

---

## Figure 4: Disease Progression and Driver Mutation Impact

**Question:** What are the key cellular and spatial changes from healthy lung through pre-invasive to invasive cancer, and how do driver mutations shape the TME?

| Panel | Content | Viz type | Data source |
|-------|---------|----------|-------------|
| **a** | Cell type abundance trajectory across disease stages (Healthy → COPD → AAH → AIS → MIA → IAC → LUAD_advanced → LUSC) | Connected dot plot / line plot with CI | healthy_lung, ggo-imc, Sorin, LUAD-vs-LUSC study |
| **b** | Immune checkpoint expression escalation by stage | Heatmap (checkpoints x cell types x stages) | ggo-imc (PD1, TIM3, CTLA4, LAG3, VISTA), public datasets |
| **c** | Spatial niche remodeling: domain composition changes across progression | Alluvial/Sankey diagram | ggo-imc UTAG, healthy_lung UTAG |
| **d** | Driver mutation impact on TME (EGFR vs KRAS vs TP53) | Forest plot of differential cell type abundance by mutation | ggo-imc mutations, Driver mutations study (Nat Commun 2025, n=157) |
| **e** | LUAD vs LUSC vs SCLC TME comparison | Side-by-side archetype prevalence + key cell type differences across histologies | LUAD-vs-LUSC study (Nat Commun 2025, n=102), Sorin, TRACERx, SCLC studies |

**Methods:** Jonckheere-Terpstra trend test for monotonic changes. Logistic regression for mutation-TME associations (multivariate).

**Stats:** BH-corrected pairwise Wilcoxon between adjacent stages. Spearman correlation of checkpoint expression with ordinal stage. Mann-Whitney for mutation comparisons.

---

## Figure 5: Clinical Translation and Spatial Biomarkers

**Question:** Can spatial TME features predict clinical outcomes and guide therapeutic strategies?

| Panel | Content | Viz type | Data source |
|-------|---------|----------|-------------|
| **a** | Survival stratification by tissue archetype | KM curves (log-rank) + Cox forest plot | Sorin (survival data, n=416), ggo-imc, TRACERx |
| **b** | Spatial biomarker model: ROC curves comparing spatial features vs cell proportions vs clinical features alone | ROC + AUROC comparison | CODEX immunotherapy cohort, Sorin |
| **c** | Top predictive spatial features (SHAP) | SHAP beeswarm plot (top 15 features) | Model from 5b |
| **d** | TLS as prognostic biomarker: density across stages + association with outcome | Violin (TLS density by stage) + KM by TLS status | ggo-imc, Sorin, CODEX |
| **e** | Therapeutic target mapping: archetype-to-therapy matching | Bubble plot (archetypes x targets, size=patient fraction, color=expression) | Integrated atlas |

**Methods:** XGBoost classifier with SHAP. lifelines for KM + Cox PH. TLS identification from B cell + Tfh co-localization in UTAG domains.

**Stats:** DeLong test for AUROC comparison. Log-rank for KM. 5-fold cross-validation. BH correction.

---

## Supplementary Figures

### S1: Data Quality and Processing
- Per-sample cell counts across cohorts (bar)
- Marker intensity distributions pre/post normalization (violin)
- UMAP pre/post Harmony (batch colored)
- Integration benchmark: kBET, LISI, silhouette across methods
- Segmentation quality examples (CellPose masks on IMC images)

### S2: Cell Type Annotation Validation
- Marker x cell type dot plots per cohort (separately)
- Automated vs manual label confusion matrix
- Panel concordance for dual-panel samples (ggo-imc G vs H)
- Subsampling stability (bootstrap CI)

### S3: T Cell Detailed Analysis
- T cell sub-cluster UMAP (8-10 subtypes) with marker overlays
- Functional marker heatmap (activation/exhaustion/memory)
- CD8/CD4 ratio across stages (violin)
- Treg spatial localization: distance to nearest tumor cell
- Cytokine+ lymphocyte density by pathology stage

### S4: Myeloid Compartment Analysis
- Myeloid sub-cluster UMAP (macrophage, monocyte, DC, neutrophil, mast, NK)
- M1/M2 markers by stage (CD163, CD206, CD68)
- Cytokine profiles by subtype (IFNg, IL1, IL12, IL17, IL23 from Panel H)
- PMN-MDSC density by stage (VISTA, CD66b)

### S5: Stromal and Vascular Analysis
- Endothelial (CD31+) density by stage
- CAF subtypes and co-localization with tumor
- Vascularization domain changes
- EMT gradient (PanCK vs Vimentin)
- Basement membrane (ColTypeIV) integrity

### S6: Spatial Domain Detailed Characterization
- Full UTAG domain taxonomy (all domains, both panels)
- Domain-level marker expression profiles
- Domain network connectivity (adjacency)
- Domain stability analysis (different radii)

### S7: Cross-Cohort Integration Details
- Per-cohort UMAP (each separately)
- Shared marker expression correlation (scatter per marker)
- UINMF partial overlap integration performance
- Imputed vs measured marker comparison

### S8: Mutation-TME Associations Extended
- Full oncoplot (EGFR, KRAS, TP53, BRAF, STK11, KEAP1)
- Per-mutation cell type density forest plots (all combinations)
- Smoking status impact on TME
- Mutation co-occurrence with TME archetypes

### S9: Healthy Lung Reference
- Healthy lung cell type UMAP and composition
- Proximal vs distal airway differences
- Healthy spatial domains (normal architecture baseline)
- COPD vs healthy comparison
- Secretory cell markers (MUC5AC, MUC5B, CC16, SCGB3A2, SFTPB/C)

### S10: Spatial Biomarker Extended Analysis
- Full feature importance ranking
- Cross-validation per fold (box plots)
- Individual biomarker KM curves (TLS density, CD8-tumor distance, M2/M1 ratio)
- Comparison to published signatures (Sorin DL model, TRACERx archetypes, ICI spatial sigs)

---

## Computational Methods Summary

| Method | Purpose | Package/Tool |
|--------|---------|-------------|
| Harmony / scVI | Cross-cohort integration on shared markers | harmonypy / scvi-tools |
| UINMF / StabMap | Partial marker overlap integration | LIGER / StabMap R |
| UTAG | Spatial domain discovery | ElementoLab/utag |
| Squidpy | Neighborhood enrichment, spatial graphs, co-localization | squidpy |
| Leiden | Cell type clustering + sub-clustering | scanpy |
| NMF / consensus clustering | Tissue archetype discovery | sklearn NMF / consensus |
| DPT / palantir | Epithelial transition pseudotime | scanpy / palantir |
| XGBoost + SHAP | Predictive modeling + interpretability | xgboost / shap |
| lifelines | KM, Cox PH, log-rank | lifelines |
| marsilea | Complex composite figures | marsilea |

---

## Implementation Phases

### Phase 0: Data Acquisition (1-2 weeks)
- [x] Download Sorin data from Zenodo (2.1 GB, record 7760826, corrected from plans wrong 13947395)
- [x] Download CODEX NSCLC from Zenodo (50 GB)
- [x] Download Desharnais LUAD-vs-LUSC from Zenodo (740 MB)
- [x] Download COVID lung IMC from Zenodo (978 MB)
- [x] Download Cords CAF from Zenodo v1 (13 GB, SCE format, Zenodo 7961844 v1 — includes SCE objects, metadata CSVs, corrected masks, documentation)
- [ ] Apply for TRACERx EGA access
- [ ] Zhu et al. — blocked, requires MTA from MD Anderson, skipping for now
- [ ] Standardize all into `adata.ingested.h5ad` per Architecture.md
- [ ] Convert Cords SCE RDS files to AnnData (anndata2ri or zellkonverter)
- [x] Unzip Sorin `LungData.zip` and Desharnais `LUSC_LUADzenodo.zip` on cayuga
- [x] Inspect Sorin/Desharnais `.mat` file structure (scipy.io.loadmat probe script)
- [x] Write MATLAB→AnnData converter for Walsh lab pipeline (covers both Sorin + Desharnais)
- [x] Unzip Cords `SingleCellExperiment_Objects.zip` on cayuga
- [x] Convert Cords SCE `.rds` → h5ad via R zellkonverter (needs 96 GB RAM, R 4.4.1 on cayuga)
- [x] Extract healthy lung subset from COVID IMC h5ad
- [ ] Assess marker panel overlap across all datasets for harmonization feasibility
- [ ] Harmonize marker names across datasets (e.g. CD8a→CD8, PanCK→Pancytokeratin)
- [ ] Move Cords spatial coords from obs to obsm['spatial']

### Phase 1: Integration & Cell Typing (2-3 weeks)
- Panel harmonization: identify shared markers, arcsinh normalize
- Integration benchmark on shared marker space (Harmony, ComBat, scVI)
- Unified cell type annotation across cohorts
- Sub-clustering of major lineages
- UTAG spatial domain computation for new datasets

### Phase 2: Analysis & Figure Generation (3-4 weeks)
- Main figures 1-5
- Statistical testing with BH correction
- Spatial archetype classification
- Predictive model training
- Disease progression analysis

### Phase 3: Supplementary & Polish (2-3 weeks)
- Supplementary figures S1-S10
- Publication-quality figure formatting (300 DPI, Helvetica/Arial 5-8pt, color-blind safe)
- Manuscript writing

---

## Key Analysis Decisions

- **Integration strategy:** Harmony on arcsinh-transformed shared markers (~12-15 core: CD3, CD4, CD8, CD68, CD163, PanCK, Ki67, aSMA, CD31, CD20, FoxP3, HLADR). UINMF/StabMap for non-overlapping markers.
- **Normalization:** `arcsinh(x/5)` for all IMC data
- **Spatial domains:** UTAG (ElementoLab) -- already computed for internal datasets, needs computation for public data
- **Archetype classification:** NMF/consensus clustering on patient-level TME composition (analogous to LME classes in lymph_dlbcl)
- **Deconvolution:** Not applicable -- IMC is single-cell resolution

---

## Key Risks & Mitigations

| Risk | Mitigation |
|------|-----------|
| Panel overlap too small for integration | Focus on ~12-15 core shared markers; UINMF for non-shared; validate with held-out markers |
| Public data not available at single-cell level | Use aggregated stats for meta-analysis panels; limit spatial analyses to datasets with coordinates |
| TRACERx access denied | Paper still viable without it; Sorin (416 pts) provides the large LUAD cohort |
| Zhu data blocked by MTA | Paper viable without it; ggo-imc covers the precursor stages (AAH/AIS/MIA/IAC) with ~40 patients |
| GGO sample size small for progression | Combine ggo-imc + ggo_human (~50 total); Zhu unavailable (MTA); focus on internal data strength |
| SCLC data limited (mostly mouse) | Include as comparison panel; combine mouse SCLC IMC (221 ROIs) + human SCLC spatial proteomics (2025); clearly label mouse vs human |
| Different segmentation pipelines across datasets | Benchmark effect; normalize at expression level (arcsinh) |

### Blockers

None currently active. (cayuga connectivity resolved 2026-03-15)

---

## Verification

- Integration quality: kBET batch score > 0.8, LISI batch mixing > 1.5
- Cell type concordance: ARI > 0.7 across cohorts for same disease stage
- Spatial archetype reproducibility: cross-validation accuracy > 75% across cohorts
- All stats: BH-corrected p-values, significance bars per skills.md standards
- Figures: 300 DPI, color-blind safe (Okabe-Ito), lowercase panel labels (Nature style)
- Run `make lint` before any commit

---

## Reusable Templates from Sibling Projects

| What | Where | Use for |
|------|-------|---------|
| Figure infrastructure (colors, rcParams) | `projects/imc/lymph_dlbcl/scripts/figure_config.py` | Consistent figure styling |
| Archetype/LME classification | `projects/imc/lymph_dlbcl/scripts/fig2_lme_classes.py` | TME archetype discovery (Fig 3d) |
| IMC loader | `sc_tools.ingest.loaders.load_imc_sample()` | Standardize public data into h5ad |
| Panel mapper | `sc_tools.ingest.IMCPanelMapper` | Resolve marker names across panels |
| Integration benchmark | `sc_tools.bm.integration` | Cross-cohort integration evaluation |
| UTAG spatial domains | `projects/imc/ggo-imc/` | Domain computation workflow |

---

## Critical Files

| File | Role |
|------|------|
| `projects/imc/imc_luad_atlas/` | Target project directory (currently empty) |
| `projects/imc/ggo-imc/` | Primary internal data (GGO, 2 panels, cell-typed) |
| `projects/imc/healthy_lung/` | Healthy reference (UTAG domains, 28 markers) |
| `projects/imc/ggo_human/` | Legacy GGO data (cell-typed, MCD files) |
| `projects/imc/lymph_dlbcl/scripts/figure_config.py` | Figure infrastructure template (colors, rcParams) |
| `projects/imc/lymph_dlbcl/scripts/fig2_lme_classes.py` | Archetype classification template |
| `sc_tools/pipeline.py` | Phase DAG to extend |
| `sc_tools/bm/` | Integration benchmarking |
