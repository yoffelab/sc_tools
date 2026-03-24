# Phase 6: Scientific Gaps - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md -- this log preserves the alternatives considered.

**Date:** 2026-03-24
**Phase:** 06-scientific-gaps
**Areas discussed:** Pseudobulk DE design, Marker validation report, Subject-level metadata model, Panel-aware cell typing

---

## Pseudobulk DE Design Formula

| Option | Description | Selected |
|--------|-------------|----------|
| Auto-infer with override | Auto-build ~ condition + batch. User override with --formula. | ✓ |
| Explicit formula only | User must always provide formula. Less agent-friendly. | |
| You decide | Claude's discretion. | |

**User's choice:** Auto-infer with override
**Notes:** Min thresholds: >=3 subjects/group, >=10 cells/subject+celltype.

---

## DE Results Output Format

| Option | Description | Selected |
|--------|-------------|----------|
| Per-celltype CSV + summary JSON | One CSV per type. Standard for downstream tools. | ✓ |
| Single wide CSV | All types in one file. | |
| AnnData varm layer | Store in adata. Non-standard for DE. | |

**User's choice:** Per-celltype CSV + summary JSON
**Notes:** Columns: gene, log2FC, pvalue, padj, baseMean.

---

## Marker Validation Report

| Option | Description | Selected |
|--------|-------------|----------|
| HTML report with dotplot | Extend existing report system. Top 5 markers/type. Flag low expression. | ✓ |
| Static PNG figures | Not integrated with report system. | |
| You decide | Claude's discretion. | |

**User's choice:** HTML report with dotplot
**Notes:** Flagging threshold default 0.1. Informational only.

---

## Subject-Level Metadata

| Option | Description | Selected |
|--------|-------------|----------|
| Warn + validate, don't block | Warn single-sample. Enforce multi-sample. Check confounding. | ✓ |
| Always require subject_id | Breaks single-sample workflows. | |
| Add at registration only | Misses ingestion-time checks. | |

**User's choice:** Warn + validate, don't block
**Notes:** subject_id is obs column convention, not schema migration.

---

## Panel-Aware Cell Typing

| Option | Description | Selected |
|--------|-------------|----------|
| Warn and restrict, allow override | Auto-restrict to sctype/custom_gates. --force-method override. | ✓ |
| Block whole-transcriptome | Hard block. No override. | |
| Warn only | Don't restrict. Agents might ignore. | |

**User's choice:** Warn and restrict, allow override
**Notes:** n_vars < 1000 threshold. Decision logged in provenance.

---

## Claude's Discretion

- PyDESeq2 wrapper details
- Dotplot rendering library choice
- Confounding detection algorithm
- Panel dispatch implementation (guard clause vs config)
- DE command group placement
- Test fixture design

## Deferred Ideas

- Cell-typing benchmarking (v2 ADV-06)
- Trajectory / RNA velocity (v2 ADV-03)
- Bootstrap uncertainty for DE
- Volcano plots for DE results
