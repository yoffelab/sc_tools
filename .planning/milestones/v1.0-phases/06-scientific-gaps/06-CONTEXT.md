# Phase 6: Scientific Gaps - Context

**Gathered:** 2026-03-24
**Status:** Ready for planning

<domain>
## Phase Boundary

sc_tools supports pseudobulk DE, marker validation, subject-level metadata, and panel-aware cell typing for rigorous multi-sample analysis. Requirements SCI-01 through SCI-04. This is the first phase adding new scientific analysis capabilities (not CLI infrastructure). Creates `sc_tools.tl.de` module, expands marker validation reporting, adds subject_id enforcement at ingestion, and adds panel detection to cell typing dispatch.

</domain>

<decisions>
## Implementation Decisions

### Pseudobulk DE (SCI-01)
- **D-01:** Auto-infer design formula with override. Auto-build `~ condition + batch` from obs metadata. Detect batch covariates from `library_id`/`batch` columns. User can override with `--formula` flag for custom covariates (e.g., `~ condition + sex + age`).
- **D-02:** Aggregation by `subject_id + celltype`. Minimum thresholds: >=3 subjects per condition group, >=10 cells per subject+celltype combination for aggregation. Combinations below threshold are excluded with a warning.
- **D-03:** Per-celltype CSV output + summary JSON. One CSV per cell type with columns: `gene`, `log2FC`, `pvalue`, `padj`, `baseMean`. Output directory: `{project_dir}/results/de/`. CLIResult.data has summary stats, CLIResult.artifacts lists all CSV paths.
- **D-04:** PyDESeq2 as the DE engine. No alternatives — this is the standard for pseudobulk.

### Marker validation report (SCI-02)
- **D-05:** HTML report extending existing report system (Phase 3 templates in `sc_tools/assets/`). Dotplot of top 5 markers per assigned cell type. Flag table for types with low canonical marker expression (mean expression below configurable threshold, default 0.1).
- **D-06:** Integrated into `sct report post-celltyping`. Summary includes: n_types tested, n_flagged, total cells. Canonical marker source is from existing gene set signatures (`sc_tools/tl/gene_sets.py`) or user-provided marker file.
- **D-07:** Flagging is informational — does not fail the command. Agent decides whether to act on flags.

### Subject-level metadata model (SCI-03)
- **D-08:** Warn + validate, don't block. Single-sample projects: warn if `subject_id` missing but don't fail. Multi-sample projects (detected via `--multi-sample` flag or multiple unique `library_id` values): require `subject_id`, validate it's distinct from `library_id`.
- **D-09:** Batch-condition confounding check at registration/QC. When `subject_id` and a condition column exist, check if batch perfectly confounds condition. Warn if so — don't block.
- **D-10:** `subject_id` is an obs column convention, not a schema migration. Enforced via validation functions, not DB constraints.

### Panel-aware cell typing (SCI-04)
- **D-11:** When `n_vars < 1000`, auto-restrict to panel-validated methods: `sctype` (with custom markers) and `custom_gates`. Warn if user requests whole-transcriptome models (celltypist, scgpt, geneformer, scarches, singler).
- **D-12:** `--force-method` flag allows override with a warning logged in provenance. The warning clearly states results may be unreliable.
- **D-13:** Panel detection logged in CLIResult provenance. Decision recorded: `{panel_detected: true, n_vars: 300, restricted_methods: [...]}`.

### Claude's Discretion
- PyDESeq2 wrapper implementation details (anndata to pandas conversion, count matrix extraction)
- Dotplot rendering library (scanpy.pl.dotplot, custom matplotlib, or Plotly for HTML)
- Exact confounding detection algorithm (perfect confounding vs partial)
- How `annotate_celltypes()` dispatch changes for panel mode (decorator, guard clause, or config)
- Whether `sct de run` is a new command group or subcommand under existing group
- Test fixture design for pseudobulk (synthetic multi-subject data)

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Cell typing dispatch (modify for SCI-04)
- `sc_tools/tl/celltype/annotate.py` -- `annotate_celltypes()` main dispatch function
- `sc_tools/tl/celltype/apply.py` -- `apply_celltype_map()`
- `sc_tools/tl/celltype/_base.py` -- Base class for cell typing methods
- `sc_tools/tl/celltype/_sctype.py` -- sctype implementation (panel-compatible)
- `sc_tools/tl/celltype/_celltypist.py` -- celltypist (whole-transcriptome)
- `sc_tools/tl/celltype/__init__.py` -- Method registry

### Gene sets and markers (SCI-02 marker source)
- `sc_tools/tl/gene_sets.py` -- Gene set definitions
- `sc_tools/tl/score_signature.py` -- Signature scoring

### QC and reporting (SCI-02 report, SCI-03 validation)
- `sc_tools/qc/sample_qc.py` -- QC metrics (add subject_id checks)
- `sc_tools/qc/report.py` -- HTML report generation (extend for post-celltyping)
- `sc_tools/assets/` -- Report HTML templates
- `sc_tools/pl/heatmaps.py` -- Existing heatmap plotting

### Ingestion (SCI-03 enforcement)
- `sc_tools/ingest/loaders.py` -- Data loading entry points (add subject_id warning)

### CLI (new commands)
- `sc_tools/cli/__init__.py` -- Command registration
- `sc_tools/cli/qc.py` -- QC commands (add subject_id validation)
- `sc_tools/models/result.py` -- CLIResult, ProvenanceRecord

### Requirements
- `.planning/REQUIREMENTS.md` -- SCI-01 through SCI-04

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `annotate_celltypes()` in `tl/celltype/annotate.py` -- cell typing dispatch, needs panel guard
- `sc_tools/qc/report.py` + `sc_tools/assets/` templates -- HTML report system for marker validation report
- `sc_tools/tl/gene_sets.py` -- canonical marker definitions
- `cli_handler` + provenance system -- auto sidecar writing for new DE command
- `sc_tools/pl/heatmaps.py` -- existing plotting for dotplot/heatmap generation

### Established Patterns
- CLI commands use `cli_handler` decorator + CLIResult envelope
- HTML reports use Jinja2 templates in `sc_tools/assets/`
- Lazy imports for heavy deps (PyDESeq2 would be optional)
- `_check_deps()` for fast-fail on missing optional deps

### Integration Points
- New `sc_tools/tl/de.py` module for pseudobulk DE
- New `sct de run` CLI command (or `sct celltype de`)
- `annotate_celltypes()` gains panel detection guard
- `sct report post-celltyping` gains marker validation section
- QC functions gain `subject_id` validation checks

</code_context>

<specifics>
## Specific Ideas

- PyDESeq2 is an optional dependency — use `_check_deps(['pydeseq2'])` fast-fail pattern from Phase 3
- DE results as per-celltype CSVs are standard for downstream tools (EnrichR, fgsea, GSEA)
- Batch-condition confounding check: if all subjects in condition A are from batch 1 and all in condition B from batch 2, warn clearly
- Panel detection threshold of n_vars < 1000 covers all targeted panels (CosMx SMI ~960, IMC ~40) while excluding whole-transcriptome (Visium ~18000, scRNA ~33000)

</specifics>

<deferred>
## Deferred Ideas

- Cell-typing benchmarking (ADV-06 in v2) — compare sctype, celltypist, scArches systematically
- Trajectory / RNA velocity (ADV-03 in v2)
- Bootstrap uncertainty for DE results
- Volcano plot generation for DE results — could be a future report enhancement

</deferred>

---

*Phase: 06-scientific-gaps*
*Context gathered: 2026-03-24*
