---
type: marker
tier: knowledge
tags:
  - knowledge/markers
---

# Immune Human Protein Panel (IMC)

40-marker IMC panel for human immune profiling. Bundled in `sc_tools/data/immune_human_protein.json`.

## Panel Categories

- **T cells:** CD3, CD4, CD8a, FoxP3, PD-1, Ki67
- **B cells:** CD20, PAX5
- **Myeloid:** CD68, CD163, CD14, CD11c, HLA-DR
- **NK:** CD56
- **Structural:** PanCK, Vimentin, SMA, Collagen I
- **Vascular:** CD31
- **Nuclear:** DNA1 (Ir191), DNA2 (Ir193)

## IMC-specific Notes

- Channel resolution: `IMCPanelMapper` resolves protein name, isotope tag, or partial match
- RGB composite default: PanCK=R, CD3=G, DNA1=B
- Normalization: arcsinh(x/5) not log1p for IMC protein data

## Related

- [[imc_integration]] for integration benchmark
- [[hallmark]] for gene-level signatures
