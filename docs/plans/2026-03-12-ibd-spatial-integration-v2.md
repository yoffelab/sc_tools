---
status: active
created: 2026-03-12
project: ibd_spatial
summary: "IBD Spatial v2: M3 full 44-sample integration (3 strategies) + M4 gene imputation (kNN + Tangram)"
---

# Plan: IBD Spatial Integration v2 — M3 Full Integration + M4 Gene Imputation

## Context

v1 (M0-M2) established that cross-platform integration of CosMx and Xenium spatial
transcriptomics data is feasible. scVI achieved batch_score 0.992 on the high-plex
pair (M2, 2552 genes), Harmony led on the low-plex pair (M1, 119 genes). All formal
success criteria passed for scVI on M2.

v2 addresses two questions:
1. **M3**: Can we unify ALL 44 samples across 6 panels into one atlas? Three competing
   strategies will be compared (single embedding, hierarchical anchor, reference-based).
2. **M4**: Can we impute missing genes in low-plex panels using the scRNA-seq reference?
   Two methods compared (kNN transfer, Tangram) with held-out gene validation.

Deliverables: integrated atlas + benchmark report + publication figures.

---

## Data Summary

| Panel | N | Genes | Patients | Disease |
|-------|---|-------|----------|---------|
| CosMx 1k | 16 | 950 | 16 (CD+UC) | Ileum+Rectum |
| Xenium MT | 16 | 377 | 16 (same) | Ileum+Rectum |
| CosMx 6k | 4 | 6,175 | 4 (UC+Healthy) | Rectum |
| Xenium 5K | 4 | 5,001 | 4 (same) | Rectum |
| Xenium withseg | 4 | 377 | 4 (same) | Rectum |
| Xenium colon | 4 | 322 | 4 (same) | Rectum |

- **44 samples** total (excluding 4 noseg)
- **4 anchor patients** (I0400, I0355, I0380, I0387) measured on 5+ panels each
- **Reference**: SAHA_IBD_RNA.h5ad (~11GB, ~18k genes, ct_major_new/ct_minor_new)
- All 52 h5ad files already converted and on cayuga

---

## Step 0: Gene Intersection Analysis

**Script**: `run_m3_gene_analysis.py` | CPU, 5 min

Load var_names from all 44 samples (no X needed). Compute:
- Full 6-panel intersection (expected ~100-300 genes)
- 5-panel intersection excluding xenium_colon (322 genes may bottleneck)
- Pairwise panel intersection matrix (6x6)
- IBD canonical marker coverage per intersection

**Output**: `results/m3_gene_analysis/` — gene lists, coverage matrix, pairwise CSV

**Gate**: If 6-panel intersection < 80 genes, drop xenium_colon from Strategy 1
(those 4 patients are still covered by Xenium 5K + withseg).

---

## M3: Three Integration Strategies (run in parallel)

### Strategy 1 — Single Embedding

**Script**: `run_m3_strategy1_single.py` | GPU A40, ~2h

All 44 samples subset to shared genes (~100-300). Concatenate, QC filter, normalize.
Run 5 methods: PCA, Harmony, scVI (n_latent=8, n_hidden=32), scANVI, BBKNN.
batch_key = "panel" (6 levels, not just 2 platforms).

Reduced scVI capacity (8/32) is deliberate — prevents overfitting on ~100-300 genes.
For scANVI, use celltype_broad labels where available, "Unknown" for withseg/colon.

**Output**: `results/m3_strategy1/adata.m3_s1.h5ad`, benchmark CSV, UMAPs

### Strategy 2 — Hierarchical Anchor

**Script**: `run_m3_strategy2_hierarchical.py` | GPU A40, ~1.5h

Keep M1 (119 genes, 32 samples) and M2 (2552 genes, 8 samples) as separate
integrated tiers. Add withseg + colon to Tier A (re-run Harmony). Bridge tiers
using the 4 anchor patients who appear in both.

Two bridging sub-strategies:
- **(a) Procrustes**: Align Tier B embedding to Tier A using anchor patient centroids
- **(b) Harmony on combined**: Concatenate both tiers (subset to ~119 shared genes),
  run Harmony with batch_key="tier"

Advantage: Tier B retains 2,552 genes for downstream biology (signatures, DE, LigRec).

**Output**: `results/m3_strategy2/` — per-tier adatas + combined embedding

### Strategy 3 — Reference-Based Transfer (scANVI)

**Script**: `run_m3_strategy3_reference.py` | GPU A100 80GB, ~4h

Train scVI→scANVI on reference (SAHA_IBD_RNA.h5ad, subsample to ~100K cells).
For each spatial panel, use `SCANVI.load_query_data()` to project spatial cells
into the reference latent space using that panel's gene set (no intersection needed).

Each panel maps with its OWN genes against the reference — avoids gene intersection
bottleneck entirely. Panels with more genes (CosMx 6k: 6,175) get better mappings.

Also produces scANVI-predicted cell types for ALL cells (including withseg/colon
which lack original annotations).

**Output**: `results/m3_strategy3/` — adata with reference embeddings, saved model,
celltype predictions

**Risk**: Low-plex panels (119-377 genes) may map poorly. Mitigate by validating
prediction accuracy on panels with known labels before trusting predictions.

### Strategy Comparison

**Script**: `run_m3_compare_strategies.py` | CPU, ~1h

Compare all three on identical metrics:

| Metric | Weight | Direction |
|--------|--------|-----------|
| batch_score (panel mixing) | 25% | Higher = better |
| platform_entropy | 15% | Higher = better |
| celltype_broad_ASW | 20% | Higher = better |
| kNN label transfer accuracy | 15% | Higher = better |
| Anchor patient concordance | 15% | Higher = better |
| Disease separation (Fisher) | 10% | Lower p = better |

**Output**: `results/m3_comparison/strategy_comparison.csv`, comparison report,
3-row UMAP grid (strategy x color-by)

---

## M3 Downstream Biology (on winning strategy)

### IBD Spatial Signatures

**Script**: `run_m3_signatures.py` | CPU, ~1h

Create `metadata/ibd_signatures.json` with 4 categories:

- **IBD_Inflammation**: TNF signaling, IL17 pathway, IFNg response, inflammasome, NF-kB
- **IBD_Fibrosis**: fibrotic core (COL1A1/FN1/ACTA2), myofibroblast, ECM remodeling (MMPs)
- **IBD_Barrier**: tight junctions (TJP1/OCLN/CLDNs), mucus (MUC2/TFF3), antimicrobial (DEFA5/LYZ)
- **IBD_Immune**: Treg, Th17, plasma/IgA, inflammatory macrophage, mature DC

Score via `sc_tools.tl.score_signature()` → obsm['signature_score']. Report gene
coverage per signature (many genes will be missing in low-plex intersection).

For Strategy 2 (hierarchical): also score Tier B with full 2,552-gene signatures.

### Disease Boundary Mapping

**Script**: `run_m3_disease_boundaries.py` | CPU, ~2h

Per sample:
1. Compute local inflammation score (kNN average of IBD_Inflammation signature, k=20)
2. Threshold into inflamed vs non-inflamed spatial domains
3. Compare with sample-level disease_state labels
4. For matched patient pairs (same patient, different panels): compare boundary maps

Produces spatial inflammation maps per sample.

### Ligand-Receptor Analysis

**Script**: `run_m3_ligrec.py` | CPU, ~6h

1. Build spatial neighbor graph: `sq.gr.spatial_neighbors(n_neighs=6)`
2. Run `sq.gr.ligrec(cluster_key="celltype_broad", n_perms=1000)` per disease group
3. FDR correct (BH), extract significant pairs
4. Compare CD-specific vs UC-specific vs shared inflammatory interactions
5. For panels without celltype labels: use scANVI predictions from Strategy 3

Run per disease group (CD, UC, Healthy), not per sample, to keep computational.

### Neighborhood Enrichment

**Script**: `run_m3_nhood.py` | CPU, ~2h

1. `sq.gr.nhood_enrichment(cluster_key="celltype_broad")` per disease group
2. `sq.gr.co_occurrence(cluster_key="celltype_broad")` for co-localization
3. Per-sample Moran's I on key signatures via `sq.gr.spatial_autocorr`
4. Compare cell type co-localization patterns: UC vs Healthy vs CD

---

## M4: Gene Imputation (depends on M3 Strategy 3 reference model)

### Strategy A — kNN from scRNA-seq Reference

**Script**: `run_m4_knn_impute.py` | GPU A40, ~1h

Use M3 Strategy 3 scANVI embeddings. For each spatial cell, find k=20 nearest
reference cells in latent space, compute distance-weighted average expression for
genes NOT in the spatial panel.

**Held-out validation**: For each panel, hold out 20% of genes that ARE measured.
Run imputation without them, compare to actual:
- Per-gene Pearson/Spearman correlation
- Per-cell correlation across held-out genes
- 5 random splits, report mean +/- std

### Strategy B — Tangram

**Script**: `run_m4_tangram.py` + `run_m4_tangram_array.sh` (--array=0-43) | GPU A40, ~10min/sample

Per-sample Tangram mapping: reference cells → spatial locations, then project ALL
reference genes to spatial. Same held-out validation as kNN.

44 samples parallelized via SLURM array (~8h total GPU, ~20min wall with array).

### Imputation Comparison

**Script**: `run_m4_compare_imputation.py` | CPU, ~1h

Compare kNN vs Tangram:
- Per-gene correlation (held-out)
- Spatial autocorrelation of imputed genes (Moran's I — should be positive for
  tissue-specific genes)
- Cell type marker recovery: do imputed markers correctly identify cell types?
- Stratified by panel (low-plex panels should benefit most)

---

## Execution DAG and GPU Budget

```
Step 0: Gene analysis (CPU, 5m)
  |
  +-- S1: Single embedding (A40, 2h)  ---|
  +-- S2: Hierarchical (A40, 1.5h)   ---|-- parallel
  +-- S3: Reference-based (A100, 4h)  ---|
       |
       Compare strategies (CPU, 1h)
       |
       +-- Signatures (CPU, 1h)  ------|
       +-- Boundaries (CPU, 2h)  ------|-- parallel
       +-- LigRec (CPU, 6h)     ------|
       +-- Nhood (CPU, 2h)      ------|
       |
       +-- M4 kNN (A40, 1h)   --------|-- parallel
       +-- M4 Tangram (A40 array, 8h) -|
            |
            M4 Compare (CPU, 1h)
            |
            Report + figures (CPU, 2h)
```

**Total GPU**: ~17h (S1:2 + S2:1.5 + S3:4 + M4-kNN:1 + M4-Tangram:8)
**Total CPU**: ~16h
**Wall time with parallelization**: ~16h

---

## Scripts to Create (13 Python + 13 SBATCH + 1 metadata)

| # | Script | GPU | Time |
|---|--------|-----|------|
| 1 | `run_m3_gene_analysis.py` + `.sh` | No | 5m |
| 2 | `run_m3_strategy1_single.py` + `.sh` | A40 | 2h |
| 3 | `run_m3_strategy2_hierarchical.py` + `.sh` | A40 | 1.5h |
| 4 | `run_m3_strategy3_reference.py` + `.sh` | A100 | 4h |
| 5 | `run_m3_compare_strategies.py` + `.sh` | No | 1h |
| 6 | `run_m3_signatures.py` + `.sh` | No | 1h |
| 7 | `run_m3_disease_boundaries.py` + `.sh` | No | 2h |
| 8 | `run_m3_ligrec.py` + `.sh` | No | 6h |
| 9 | `run_m3_nhood.py` + `.sh` | No | 2h |
| 10 | `run_m4_knn_impute.py` + `.sh` | A40 | 1h |
| 11 | `run_m4_tangram.py` + `_array.sh` | A40 | 8h |
| 12 | `run_m4_compare_imputation.py` + `.sh` | No | 1h |
| 13 | `generate_m3m4_report.py` + `.sh` | No | 2h |
| - | `metadata/ibd_signatures.json` | - | - |

---

## Publication Figures

| Figure | Content | Dir |
|--------|---------|-----|
| Fig 1 | Study design: panel diagram, patient matching, gene overlap heatmap | manuscript/ |
| Fig 2 | Integration progression: M0→M3 UMAP grid, batch_score trajectory | manuscript/ |
| Fig 3 | Strategy comparison: 3-panel UMAP (S1/S2/S3) x celltype + disease + platform | manuscript/ |
| Fig 4 | IBD spatial signatures: inflammation/fibrosis/barrier score maps, CD vs UC | manuscript/ |
| Fig 5 | LigRec: significant L-R dotplot, disease-specific interactions | manuscript/ |
| Fig 6 | Imputation: held-out validation, spatial autocorrelation of imputed genes | manuscript/ |
| Supp 1 | Per-panel QC: violins of counts, genes, cell types | supplementary/ |
| Supp 2 | Full benchmark: scib metrics table, per-method UMAPs all strategies | supplementary/ |
| Supp 3 | Neighborhood enrichment heatmaps per disease | supplementary/ |
| Supp 4 | Imputation per-gene scatters (top/bottom genes) | supplementary/ |

Standards: 300 DPI, Helvetica 5-8pt, Okabe-Ito palette, scale bars, FDR p-values.

---

## Risk Mitigation

| Risk | Mitigation |
|------|------------|
| Gene intersection < 80 | Drop xenium_colon from S1 (4 patients still covered by 5K+withseg) |
| S3 reference mapping poor for low-plex | Validate on panels with known labels first; two-stage bridge if needed |
| A100 OOM on reference training | Subsample reference to 100K cells, reduce batch_size |
| Tangram fails on small samples | Skip samples < 500 cells, note in report |
| S2 anchor bridging poor | Per-celltype anchoring; if still poor, report as informational |
| LigRec too slow (44 samples x 1000 perms) | Run per disease group (3 groups), subsample to 50K/group |
| Missing celltype labels (withseg/colon) | Use S3 scANVI predictions; flag as predicted |

---

## Verification

After all steps complete:
1. All 3 strategy adatas load and contain expected obsm keys
2. Benchmark CSVs have metrics for all methods
3. IBD signatures scored with >50% gene coverage for at least 3/4 categories
4. LigRec produces significant pairs (FDR < 0.05) for at least 2 disease groups
5. Imputation validation: median held-out gene Pearson r > 0.3 for best method
6. Report HTML loads with all embedded images
7. Publication figures pass `figure_quality_check.py` hook at manuscript/ level
8. Mission.md and Journal.md updated with M3/M4 results
