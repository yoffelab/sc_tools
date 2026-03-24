# Phase 8: Multi-Omic Assembly - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-03-24
**Phase:** 08-multi-omic-assembly
**Areas discussed:** Metadata join strategy, MuData structure, Joint embedding approach, Cross-modal query API

---

## Area Selection

| Option | Description | Selected |
|--------|-------------|----------|
| Metadata join strategy (Recommended) | How should subject_id joins work across modalities? | ✓ |
| MuData structure | How should the assembled MuData be organized? | ✓ |
| Joint embedding approach | Which integration method(s) to support? | ✓ |
| Cross-modal query API | What queries should sct support? | ✓ |

**User's choice:** All four areas
**Notes:** None

---

## Metadata Join Strategy

### Q1: Incomplete modality coverage

| Option | Description | Selected |
|--------|-------------|----------|
| Outer join (Recommended) | Include all patients even with missing modalities | ✓ |
| Inner join | Only patients present in ALL modalities | |
| Configurable | Flag to choose inner/outer at assembly time | |

**User's choice:** Outer join
**Notes:** None

### Q2: Multiple samples per modality per patient

| Option | Description | Selected |
|--------|-------------|----------|
| Patient-level summary (Recommended) | Aggregate to one row per patient | |
| Sample-level rows | One row per sample across all modalities | |
| Hierarchical | Patient-level table plus per-modality sample tables | ✓ |

**User's choice:** Hierarchical
**Notes:** None

---

## MuData Structure

### Q1: Per-modality organization in MuData

| Option | Description | Selected |
|--------|-------------|----------|
| Modality-keyed (Recommended) | Keys by modality type: mdata['rna'], mdata['imc'] | |
| Sample-keyed | Keys by sample ID | |
| Platform-keyed | Keys by platform+batch | |

**User's choice:** Other — Multilayer analysis with modality-keyed base plus a wrapper that stacks and aggregates samples on different abstraction layers
**Notes:** "baseline is ['rna'], ['imc'] on a cell level. we will stack multiple mdata on patient level, and on sample level (normal within patient, tumor within patient, per slide)"

### Q2: Confirm multilayer vision

| Option | Description | Selected |
|--------|-------------|----------|
| Yes, that's right | Modality-keyed MuData with aggregation wrapper on top | ✓ |
| Not quite | Let me clarify further | |

**User's choice:** Yes — confirmed. Added: "by design the stacked mudata might require more than 3 level (ie organ level as well)"
**Notes:** N-level abstraction stacking, not hardcoded to 3 levels

### Q3: Feature space handling

| Option | Description | Selected |
|--------|-------------|----------|
| Keep separate (Recommended) | Each modality keeps its own var (features) | ✓ |
| Common gene subset | Intersect features across modalities | |
| Gene-protein mapping | Map protein markers to gene names | |

**User's choice:** Keep separate
**Notes:** None

### Q4: Wrapper API design

| Option | Description | Selected |
|--------|-------------|----------|
| Class wrapper (Recommended) | MultiOmicAtlas class wrapping MuData | ✓ |
| Utility functions | Standalone functions | |
| You decide | Claude chooses | |

**User's choice:** Class wrapper
**Notes:** None

---

## Joint Embedding Approach

### Q1: Which methods to support

| Option | Description | Selected |
|--------|-------------|----------|
| MOFA+ (Recommended) | Factor analysis, handles missing modalities | ✓ |
| MultiVI | Deep generative model from scvi-tools | ✓ |
| TotalVI | RNA + protein joint model from scvi-tools | ✓ |
| WNN | Weighted Nearest Neighbors | |

**User's choice:** MOFA+, MultiVI, TotalVI (three methods)
**Notes:** None

### Q2: Dispatch pattern

| Option | Description | Selected |
|--------|-------------|----------|
| Pluggable dispatch (Recommended) | Same pattern as annotate_celltypes() | ✓ |
| Separate commands | sct assemble mofa, sct assemble multivi, etc. | |
| You decide | Claude chooses | |

**User's choice:** Pluggable dispatch
**Notes:** None

---

## Cross-Modal Query API

### Q1: CLI command structure

| Option | Description | Selected |
|--------|-------------|----------|
| New sct assemble group (Recommended) | sct assemble build/embed/query | ✓ |
| Extend existing groups | Add to existing command groups | |
| You decide | Claude chooses | |

**User's choice:** New sct assemble group
**Notes:** None

### Q2: Essential v1 queries

| Option | Description | Selected |
|--------|-------------|----------|
| Cell type proportions (Recommended) | Per-patient proportions across modalities | ✓ |
| Modality coverage summary | Which patients have which modalities | |
| Cross-modal marker comparison | Compare shared markers across modalities | |
| Patient-level clinical correlation | Correlate features with clinical outcomes | |

**User's choice:** Cell type proportions only for v1
**Notes:** Other query types deferred to v2

---

## Claude's Discretion

- MultiOmicAtlas serialization format
- Internal method registry structure for joint embedding
- How sct assemble build discovers per-modality files
- How hierarchical metadata is stored in MuData
- Test fixture design for multi-omic scenarios
- MOFA+/MultiVI/TotalVI wrapper thickness

## Deferred Ideas

- Cross-modal marker comparison (v2 query)
- Modality coverage summary (v2 query)
- Patient-level clinical correlation (v2 query)
- WNN embedding method (v2)
