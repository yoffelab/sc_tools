---
type: method
tier: knowledge
tags:
  - knowledge/methods
---

# IMC Integration Methods

Benchmark results from ggo_human IMC project (~248K cells, 38 markers).

## Benchmark Results (2026-03-06)

| Rank | Method | Batch Score | Notes |
|------|--------|-------------|-------|
| 1 | Z-score + Harmony | 0.634 | **Selected** |
| 2 | IMC Pheno + Harmony | 0.633 | Near-identical to #1 |
| 3 | IMC Phenotyping | 0.627 | Per-ROI z-score + ComBat + BBKNN |
| 4 | scANVI | 0.621 | Requires preliminary labels |
| 5 | ComBat | 0.568 | |
| 6 | PCA (unintegrated) | 0.564 | Baseline |
| 7 | scVI | 0.555 | |
| 8 | Harmony | 0.553 | Standard Harmony on log1p |
| 9 | CytoVI | 0.486 | Specialized IMC VAE |

## Key Decision

Batch score is the primary metric for integration selection at `preprocess`. Bio metrics (ARI, NMI) are meaningful only after validated celltyping (`celltype_manual`). See [[checkpoints.gen#preprocess]] for `preprocess` requirements.

## Harmony PyTorch Fix

Direct `harmonypy.run_harmony()` call bypasses scanpy wrapper. Handles `Z_corr` as torch tensor (`.detach().cpu().numpy()`) or numpy array.

## Related

- [[scvi_visium]] for transcriptomic integration
- [[conventions]] for statistical standards
