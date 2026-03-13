---
type: hypotheses
project: lymph_dlbcl
platform: imc
---

# lymph_dlbcl Hypotheses

Scientific questions for the Lymph Node DLBCL IMC project (328 treatment-naive tumors, 52 markers, 2 panels: T2 immune + S2 stromal).

## H1: DLBCL microenvironment hierarchy

**Status:** #hypothesis/investigating
**Question:** Do five lymph node microenvironment (LME) classes (Cold, CD206 Enriched, Cytotoxic, Stromal, T cell Regulated) predict patient outcome?
**Approach:** k=30 kNN spatial community detection, three-tier LME assignment (patient_tme CSV, TME cluster mapping, k-means fallback).

## H2: B cell subtype composition

**Status:** #hypothesis/investigating
**Question:** Does B cell subtype composition differ across LME classes and predict outcome?
**Approach:** Cell typing with 52-marker panel, per-LME composition analysis.

## H3: Spatial community structure

**Status:** #hypothesis/investigating
**Question:** Do spatial communities (k=30 kNN) correlate with LME class assignment?
**Approach:** Community detection on spatial graph, mapping to pre-defined LME classes.

## H4: Vessel density by LME class

**Status:** #hypothesis/investigating
**Question:** Does vessel density (CD31+ area) vary across LME classes?
**Approach:** Vessel segmentation from CD31 channel, per-LME quantification.
**Related:** [[immune_human_protein]]
