---
type: method
tier: knowledge
tags:
  - knowledge/methods
---

# scVI for Visium Data

Best practices for scVI integration on Visium spatial transcriptomics.

## Recipe

The `sc_tools.pp.preprocess(modality="visium")` recipe uses:
1. Raw counts (no normalization before scVI)
2. Seurat v3 HVG selection (~3000 genes)
3. scVI model training with `batch_key` per library
4. Embedding stored in `obsm['X_scvi']`

## When to Use

- Multiple Visium libraries needing batch correction
- Default integration method for transcriptomic spatial data
- Prefer over Harmony when batch effects are strong

## Alternatives

For non-VAE workflows: `normalize_total` + `log1p` + Seurat HVG + Harmony on PCA.

## Related

- [[imc_integration]] for protein-level data
- [[conventions]] for embedding naming conventions
