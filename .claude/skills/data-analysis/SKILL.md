---
name: data-analysis
tier: 3
description: "[Tier 3 — Specialized] Systematic data analysis: explore first, hypothesize, validate statistically. Use when exploring a new dataset, diagnosing unexpected results, or designing an analysis pipeline for sc_tools projects."
allowed-tools: Read Bash Glob Grep
---

# Data Analysis

## Core Principle

Never jump to conclusions from raw numbers. Explore → hypothesize → validate → report.

## Phase 1: Exploration (EDA)

Before any modeling or statistical testing:

1. **Shape and types:** `adata.shape`, `adata.obs.dtypes`, `adata.var.dtypes`
2. **Missingness:** `adata.obs.isna().sum()` — understand what is absent
3. **Distributions:** histograms of `total_counts`, `n_genes_by_counts`, `pct_counts_mt`
4. **Outliers:** MAD-based outlier detection before removing anything
5. **Batch structure:** how many samples per `library_id`? per `batch`?

```python
import scanpy as sc
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], jitter=0.4)
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color="pct_counts_mt")
```

## Phase 2: Hypothesis Formation

Write down what you expect to see BEFORE running statistical tests. Prevents p-hacking.

- "Tumor samples will have higher `pct_counts_mt` than normal"
- "Gene signature X will correlate with clinical outcome Y"

## Phase 3: Statistical Validation

Follow sc_tools standards (skills.md §10):

- **FDR:** Benjamini-Hochberg always; never report uncorrected p-values
- **Group comparisons:** 2 groups → pairwise Wilcoxon / t-test; >2 groups → 1-vs-rest (default) or all-pairwise
- **Effect size:** report alongside p-value (Cohen's d, fold change, AUC)
- **Sample size:** check N per group before testing; do not test groups with N < 5
- **Multiple comparisons:** adjust across all tests run in one analysis block

```python
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

stats, pvals = zip(*[mannwhitneyu(a, b) for a, b in pairs])
_, padj, _, _ = multipletests(pvals, method="fdr_bh")
```

## Phase 4: Visualization

- Show the data, not just summary statistics
- Violin + strip plot > bar chart for distributions
- Spatial plots: always include scale bar; colorblind-safe palettes (Okabe-Ito default)
- Label significance bars with adjusted p-values: `*` < 0.05, `**` < 0.01, `***` < 0.001

## sc_tools Workflow

| Analysis type | sc_tools function | Output |
|---------------|-------------------|--------|
| QC exploration | `sc_tools.qc.calculate_qc_metrics()` | `obs` columns |
| Gene scoring | `sc_tools.tl.score_signature()` | `obsm['signature_score']` |
| Enrichment | `sc_tools.tl.run_ora()` | DataFrame with padj |
| Spatial autocorrelation | `sc_tools.qc.spatially_variable_genes()` | `var` with Moran's I |
| Deconvolution | `sc_tools.tl.deconvolution()` | `obsm['cell_type_proportions']` |

## Common Mistakes

1. **Normalizing before QC** — always QC on raw counts
2. **Testing on the same data used to discover the hypothesis** — use held-out data or pre-register
3. **Treating clusters as independent groups** — they are not independent; correct accordingly
4. **Plotting mean ± SEM without showing N** — always report sample size
5. **Using jet/rainbow colormaps** — perceptually misleading; use viridis/inferno or Okabe-Ito
6. **Subsetting adata before checking if the subset is representative** — check composition first

## Reporting

Include in every analysis summary:
- N (spots, cells, patients) per group
- Test used and its assumptions
- Effect size and direction
- Adjusted p-value with correction method
- Figure number and exact caption
