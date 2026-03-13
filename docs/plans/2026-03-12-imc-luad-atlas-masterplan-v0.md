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

| Dataset | Citation | Cancer type | N patients | Markers | Data access | Download |
|---------|----------|------------|-----------|---------|-------------|----------|
| **Sorin et al.** | Nature 2023 | LUAD (5 histological patterns) | 416 | 35+ | Public | [Zenodo 13947395](https://zenodo.org/records/13947395) — CSV + MCD |
| **Zhu et al.** | Cancer Cell 2025 | LUAD precursors (TIM-3 focus) | 123 samples, 1655 ROIs | IMC + WES + spatial tx | Check Synapse / supp | May be on Synapse; check Cancer Cell data availability |
| **Cords et al. (CAF)** | Cancer Cell 2024 | NSCLC | 1,070 | 45 | Check Synapse / supp | Largest IMC NSCLC cohort — 11 CAF phenotypes, 4 prognostic groups |
| **LUAD vs LUSC** | Nat Commun 2025 | NSCLC (balanced LUAD/LUSC) | 102 | IMC | Check supp | Histology-specific TME differences |
| **Driver mutations** | Nat Commun 2025 | LUAD | 157 | IMC | Check supp | KRAS/EGFR/TP53 impact on TME |
| **TRACERx IMC** | Cancer Discov 2024 | NSCLC (early-stage) | 81 (198 samples) | 2 panels, 2.3M cells | Restricted (EGA DAC) | Apply: EGAC00001000632; PHLEX test data public on [Zenodo 7973724](https://zenodo.org/records/7973724) |
| **CODEX NSCLC** | J Transl Med 2024 | NSCLC (immunotherapy) | ~20+ | 36 (CODEX) | Public | [Zenodo 10258578](https://zenodo.org/records/10258578) — h5ad + TIFF (41 GB) |
| **MultiModal NSCLC** | Synapse | NSCLC (ICI response) | TBD | IMC + radiology + genomics | Public (Synapse) | [syn26642505](https://www.synapse.org/Synapse:syn26642505) |

### Public — Secondary (SCLC + reference)

| Dataset | Citation | Notes |
|---------|----------|-------|
| **SCLC IMC** | BMC 2025 | 221 ROIs, 37 markers, mouse models — include for pan-lung SCLC comparison |
| **SCLC spatial proteomics** | 2025 | Human SCLC, multiregional (center, margin, peri-tumor), survival associations |
| **Zhu mouse precancers** | Adv Sci 2026 | 284 IMC ROIs, 1.4M cells, mouse LUAD progression models |
| **COVID lung IMC** | Nature 2021 | ElementoLab, 36 markers, healthy lung reference, [Zenodo 4139443](https://zenodo.org/records/4139443) |
| **ICI spatial signatures** | Nat Genet 2025 | 234 NSCLC, spatial proteomics (n=67) + spatial tx (n=131), response/resistance sigs |
| **LuCA scRNA atlas** | cellxgene | 1.2M cells, 309 patients — reference for cell type annotation transfer |
| **HTAN lung** | Synapse | Multi-modality lung cancer data; check for IMC-specific datasets |

### Data Download Strategy

All downloads go to HPC: `/athena/elementolab/scratch/juk4007/imc_atlas_data/`

**Immediate (public, no restrictions):**
```bash
# Sorin et al. — Zenodo
zenodo download --doi 10.5281/zenodo.7383627  # or wget
# CODEX NSCLC — Zenodo (41 GB)
wget -c https://zenodo.org/records/10258578/files/anndata_4301_v4_annotated.h5ad
# COVID lung IMC (ElementoLab) — Zenodo
wget -c https://zenodo.org/records/4139443/files/*
# Synapse MultiModal NSCLC
synapse get syn26642505
```

**Requires access application:**
- TRACERx: EGA DAC (ctc.tracerx@ucl.ac.uk)
- HTAN: Check individual center submissions on Synapse

**Check supplementary / contact authors:**
- Zhu et al. 2025, Cords et al. 2024, LUAD-vs-LUSC 2025, Driver mutations 2025
- Try Synapse first, then check paper Data Availability sections

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
- Download Sorin data from Zenodo to HPC (`/athena/elementolab/scratch/juk4007/`)
- Download CODEX NSCLC from Zenodo
- Apply for TRACERx EGA access
- Track down Zhu et al. 2025 data (check Cancer Cell supplementary / contact authors)
- Check data availability for CAF (Cancer Cell 2024), LUAD-vs-LUSC (Nat Commun 2025), Driver mutations (Nat Commun 2025)
- Standardize all into `adata.ingested.h5ad` per Architecture.md

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

## Key Risks & Mitigations

| Risk | Mitigation |
|------|-----------|
| Panel overlap too small for integration | Focus on ~12-15 core shared markers; UINMF for non-shared; validate with held-out markers |
| Public data not available at single-cell level | Use aggregated stats for meta-analysis panels; limit spatial analyses to datasets with coordinates |
| TRACERx access denied | Paper still viable without it; Sorin (416 pts) provides the large LUAD cohort |
| GGO sample size small for progression | Combine ggo-imc + ggo_human + Zhu precursor data to increase N |
| SCLC data limited (mostly mouse) | Include as comparison panel; combine mouse SCLC IMC (221 ROIs) + human SCLC spatial proteomics (2025); clearly label mouse vs human |
| Different segmentation pipelines across datasets | Benchmark effect; normalize at expression level (arcsinh) |

---

## Verification

- Integration quality: kBET batch score > 0.8, LISI batch mixing > 1.5
- Cell type concordance: ARI > 0.7 across cohorts for same disease stage
- Spatial archetype reproducibility: cross-validation accuracy > 75% across cohorts
- All stats: BH-corrected p-values, significance bars per skills.md standards
- Figures: 300 DPI, color-blind safe (Okabe-Ito), lowercase panel labels (Nature style)
- Run `make lint` before any commit

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
