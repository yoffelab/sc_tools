---
type: hypotheses
project: ggo_visium
platform: visium
---

# ggo_visium Hypotheses

Scientific questions for the GGO Visium project (29,952 spots across Normal, GGO, Part-Solid, Solid).

## H1: TLS heterogeneity across tumor progression

**Status:** #hypothesis/investigating
**Question:** Do tertiary lymphoid structures differ in composition and spatial organization across lung tumor progression stages (Normal, GGO, Part-Solid, Solid)?
**Approach:** Spatial neighborhood analysis of B cell / T cell clusters, deconvolution via Cell2location + Tangram. 1-vs-rest comparisons across tumor stages.
**Related:** [[hallmark]], [[scvi_visium]]

## H2: Macrophage-proliferative spatial correlation

**Status:** #hypothesis/investigating
**Question:** Do macrophage subtypes show distinct spatial localization patterns relative to proliferative tumor programs?
**Approach:** Colocalization analysis using signature scores (Moran I, Pearson correlation). Macrophage core signature + Ki67/MKI67 scoring.
**Related:** `figures/manuscript/macrophage_localization/`

## H3: Neutrophil-cytotoxic T cell colocalization

**Status:** #hypothesis/investigating
**Question:** Are SLC16A3+ neutrophils spatially colocalized with cytotoxic T cells?
**Approach:** Colocalization via spatial multipage analysis. SLC16A3 as neutrophil marker.
**Related:** `figures/manuscript/neutrophil_cytotoxic_tcell_localization/`

## H4: Process colocalization

**Status:** #hypothesis/investigating
**Question:** Which biological processes (from signature scores) show spatial correlation with each other?
**Approach:** Moran I and Pearson correlation among all signature scores in `obsm['signature_score']`.
