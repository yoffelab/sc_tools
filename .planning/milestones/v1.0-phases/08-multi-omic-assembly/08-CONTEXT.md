# Phase 8: Multi-Omic Assembly - Context

**Gathered:** 2026-03-24
**Status:** Ready for planning

<domain>
## Phase Boundary

Independently processed modalities (scRNA, IMC, Visium, Xenium) can be assembled into a unified MuData object for cross-modal patient-level analysis. Requirements MOM-01 through MOM-04. This phase creates a new `sc_tools/assembly/` module with a `MultiOmicAtlas` class wrapper around MuData, implements patient/subject metadata joins via outer join, adds pluggable joint embedding (MOFA+, MultiVI, TotalVI), and exposes a new `sct assemble` CLI command group.

</domain>

<decisions>
## Implementation Decisions

### Metadata join strategy (MOM-01)
- **D-01:** Outer join on `subject_id`. Include all patients even with missing modalities — missing modalities are null in MuData. Filtering to complete cases is easy downstream.
- **D-02:** Hierarchical metadata structure. Patient-level table plus per-modality sample tables. Supports queries at multiple granularities: cell level, sample level (normal vs tumor, per slide), patient level, and potentially organ level.
- **D-03:** N-level abstraction stacking. The stacking is not limited to 3 fixed levels — supports arbitrary abstraction layers (cell, sample, patient, organ, cohort, etc.) defined by grouping columns in obs metadata.

### MuData structure (MOM-02)
- **D-04:** Modality-keyed MuData. Keys by modality type: `mdata['rna']`, `mdata['imc']`, `mdata['visium']`, `mdata['xenium']`. Multiple samples per modality are concatenated within each slot. Standard MuData convention.
- **D-05:** Feature spaces kept separate. Each modality keeps its own var (features). Cross-modal analysis happens through shared obs (patient/sample metadata) and joint embeddings in obsm — not by merging feature spaces.
- **D-06:** Class wrapper — `MultiOmicAtlas`. Holds MuData + patient/sample metadata tables. Methods like `.patient_view(patient_id)`, `.sample_view(sample_id)`, `.celltype_proportions()`. Encapsulates the multilayer aggregation logic.

### Joint embedding approach (MOM-03)
- **D-07:** Three methods supported: MOFA+ (factor analysis, handles missing modalities), MultiVI (deep generative, RNA+protein), TotalVI (RNA+protein joint model). MOFA+ is the recommended default for multi-platform spatial data.
- **D-08:** Pluggable dispatch pattern. Same architecture as `annotate_celltypes()` — method registry, runtime selection via `--method` flag. Consistent with existing codebase. Easy to add new methods later.

### Cross-modal query API (MOM-04)
- **D-09:** New `sct assemble` CLI command group with subcommands: `sct assemble build` (build MuData from per-modality h5ad files), `sct assemble embed` (joint embedding), `sct assemble query` (cross-modal queries). Keeps multi-omic commands together.
- **D-10:** V1 query scope: cell type proportions per patient across modalities. Core use case: "Show me B cell fractions for patient X in RNA vs IMC vs Visium." Additional query types (modality coverage, marker comparison, clinical correlation) deferred to v2.

### Claude's Discretion
- How `MultiOmicAtlas` serializes to disk (MuData h5mu format, or separate h5ad + metadata JSON)
- Internal structure of the method registry for joint embedding (enum-based vs class-based dispatch)
- How `sct assemble build` discovers and validates per-modality h5ad files (directory scan vs manifest file vs explicit paths)
- How hierarchical metadata tables are stored in MuData (uns dict, separate DataFrames, or muon's built-in obs mechanism)
- Test fixture design for multi-omic scenarios (synthetic multi-modality data)
- Whether MOFA+, MultiVI, TotalVI wrappers are thin (delegate to library) or thick (handle data prep + postprocessing)

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Subject-level metadata (Phase 6 output — join key)
- `sc_tools/qc/metadata.py` — `validate_subject_metadata()` — validates subject_id presence and distinctness
- `sc_tools/tl/de.py` — Uses `subject_id` for pseudobulk aggregation — same join key for multi-omic assembly

### Existing modality handling
- `sc_tools/ingest/loaders.py` — Platform-specific loaders (Visium, Xenium, CosMx, IMC) — defines what per-modality AnnData looks like
- `sc_tools/biodata.py` — `PlatformSpec` registry — maps platforms to modalities
- `sc_tools/pipeline.py` — `STANDARD_PHASES` with `PhaseSpec` — pipeline DAG structure

### CLI infrastructure (integration points)
- `sc_tools/cli/__init__.py` — CLI app registration pattern for new command groups
- `sc_tools/cli/estimate.py` — Recent Phase 7 example of new command group
- `sc_tools/models/result.py` — CLIResult envelope for structured output

### IO Gateway (Phase 7 output — memory-safe loading)
- `sc_tools/io/gateway.py` — IOGateway with DataTier — use for loading large per-modality h5ad files
- `sc_tools/io/estimate.py` — Memory estimation — use before assembling large multi-modal datasets

### Registry (patient/subject tracking)
- `sc_tools/registry.py` §1369-1449 — `add_subject()`, `get_subject()`, `list_subjects()` — existing patient registry

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `validate_subject_metadata()`: Already validates subject_id — reuse for assembly input validation
- `IOGateway`: Use tiered loading from Phase 7 to safely load large per-modality h5ad files
- `cli_handler`: Extend for new `sct assemble` commands with dry-run and memory guard support
- Registry patient model: `add_subject()`, `list_subjects()` — patient-level metadata store
- `annotate_celltypes()` dispatch pattern: Model for pluggable joint embedding dispatch

### Established Patterns
- Pluggable method dispatch via class registry (celltype annotation in `sc_tools/tl/celltype/`)
- CLI command groups registered in `sc_tools/cli/__init__.py`
- CLIResult envelope for all command output
- Provenance sidecars for all data-producing commands

### Integration Points
- New `sc_tools/assembly/` module for MultiOmicAtlas and joint embedding
- New `sc_tools/cli/assemble.py` for `sct assemble` command group
- Register `assemble_app` in `sc_tools/cli/__init__.py`
- MuData (muon package) as new dependency in pyproject.toml

</code_context>

<specifics>
## Specific Ideas

- User wants multilayer analysis with N-level abstraction stacking (cell → sample → patient → organ → cohort), not hardcoded to 3 levels
- Sample-level grouping includes biological distinctions within patients (normal vs tumor tissue, per slide)
- MOFA+ is preferred default because it handles missing modalities well (matches outer join strategy)
- Cell type proportions across modalities per patient is the core v1 query

</specifics>

<deferred>
## Deferred Ideas

- Cross-modal marker comparison (e.g., CD20 protein in IMC vs MS4A1 gene in RNA) — v2 query type
- Modality coverage summary table — v2 query type
- Patient-level clinical correlation (survival, treatment response) — v2 query type, builds on Phase 6 pseudobulk DE
- WNN (Weighted Nearest Neighbors) as additional embedding method — v2

</deferred>

---

*Phase: 08-multi-omic-assembly*
*Context gathered: 2026-03-24*
