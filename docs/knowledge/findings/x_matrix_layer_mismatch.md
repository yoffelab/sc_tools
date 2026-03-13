---
type: finding
tier: knowledge
tags:
  - finding/validated
  - knowledge/findings
created: 2026-03-07
projects:
  - lymph_dlbcl
validated_in:
  - lymph_dlbcl
---

# X Matrix / Layer Mismatch in Seurat-converted AnnData

## Summary

When converting Seurat objects to AnnData (e.g. via SeuratDisk), the `X` matrix may contain all zeros while actual expression data resides in `layers['raw']` or `layers['data']`. This is a silent data corruption that passes shape validation but produces meaningless downstream results.

## Evidence

Discovered in lymph_dlbcl p4 checkpoint (2026-03-07). All 47 Seurat h5ad objects had `X = 0` with expression in `layers['raw']`. Required fixing 8 analysis scripts to use `layers['raw']` fallback.

## Recommendation

Always check `adata.X.sum()` after loading Seurat-converted h5ad files. If zero, check `adata.layers` for the actual data. The `sc_tools.validate` module should be extended to catch this pattern.

## Related

- [[checkpoints.gen]] for required X format per phase
- [[conventions]] for validation standards
