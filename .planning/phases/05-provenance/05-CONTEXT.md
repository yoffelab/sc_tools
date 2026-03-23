# Phase 5: Provenance - Context

**Gathered:** 2026-03-23
**Status:** Ready for planning

<domain>
## Phase Boundary

Every CLI output has a traceable lineage back to its inputs, parameters, and environment. Requirements PRV-01 through PRV-05. This phase expands the minimal `Provenance` model from Phase 2, adds automatic sidecar writing to `cli_handler`, creates two new commands (`sct provenance show`, `sct provenance trace`), embeds provenance in `adata.uns` for portability, and threads `random_state` through clustering functions for Leiden reproducibility.

</domain>

<decisions>
## Implementation Decisions

### Sidecar writing strategy
- **D-01:** Auto via `cli_handler`. After successful execution, if `result.artifacts` is non-empty, write a `.provenance.json` sidecar next to each artifact automatically. No opt-out flag — provenance is always written.
- **D-02:** Sidecar naming: `{original_filename}.provenance.json` (e.g., `adata.normalized.h5ad.provenance.json`).
- **D-03:** Only written on success — error results don't generate sidecars.

### Provenance storage (hybrid)
- **D-04:** Hybrid approach: always write sidecar `.provenance.json` files (works for all output types — h5ad, HTML, CSV). Additionally, when the output artifact is an h5ad file, embed provenance in `adata.uns['sct_provenance']` so it travels with the data and survives file moves.
- **D-05:** Sidecars are the canonical provenance chain that `sct provenance trace` follows. `adata.uns` is the portable backup — trace falls back to it when a sidecar is missing.

### Provenance model fields
- **D-06:** Expand the existing `Provenance` Pydantic model in `models/result.py` with full fields: `command`, `params` (all flags as-passed including defaults), `inputs` (list of input file records), `sc_tools_version`, `timestamp`, `runtime_s`, `peak_memory_mb`.
- **D-07:** All CLI flags serialized as-passed — includes defaults. Agent can replay the exact command from provenance.
- **D-08:** Input file records include: `path` (relative to project root), `path_type` ("relative"), `sha256` (computed at read time), `size_bytes`.

### File portability
- **D-09:** Store paths relative to project root in sidecars (using `--project-dir` from CLI context or `.` default).
- **D-10:** When `sct provenance trace` encounters a missing input file path: (1) check `adata.uns['sct_provenance']` for embedded provenance, (2) search the project directory for a file with matching SHA256, (3) if found at a new location, continue trace and note the relocation. If not found, mark as "origin (no provenance)" and stop that branch.

### Lineage trace design
- **D-11:** `sct provenance trace <file>` follows sidecar input references recursively. Each sidecar's `inputs[].path` points to the upstream file — find its sidecar, read its inputs, recurse. Stop when no sidecar/uns provenance exists (raw data origin).
- **D-12:** Output is a flat list of lineage steps ordered chronologically (oldest first). Each step includes: file path, command (or null for origin), timestamp, and a note field for special cases (relocated, missing, origin).
- **D-13:** `sct provenance show <file>` simply reads and displays the sidecar (or `adata.uns` fallback) — no recursion.

### Leiden reproducibility (PRV-05)
- **D-14:** Thread `random_state` parameter through all functions that call `sc.tl.leiden()`. Add `random_state=0` as default parameter to `compute_integration_metrics` and any other function in the call path.
- **D-15:** Record `resolution` and `random_state` in provenance sidecar params. Test: identical inputs + identical params = identical cluster labels.

### Claude's Discretion
- How to compute SHA256 efficiently for large h5ad files (full file vs header-only)
- Internal structure of the expanded `Provenance` model (nested Pydantic models vs flat dict)
- How `cli_handler` passes input file info to the provenance (decorator metadata, return value augmentation, or explicit `inputs` field on CLIResult)
- How `adata.uns['sct_provenance']` is structured (full provenance dict or subset)
- Depth limit for trace recursion (or unlimited with cycle detection)
- How peak_memory_mb is measured (tracemalloc, resource module, or platform-specific)

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Provenance model (Phase 2 output — expand this)
- `sc_tools/models/result.py` -- `Provenance` class (minimal: command, timestamp, version) and `CLIResult` with `artifacts` field
- `sc_tools/models/__init__.py` -- model exports

### CLI handler (injection point)
- `sc_tools/cli/__init__.py` -- `cli_handler` decorator, `_emit()`, `_render_rich()`, `_state` dict
- `sc_tools/cli/qc.py` -- example command populating `artifacts` (report_path)
- `sc_tools/cli/preprocess.py` -- example command with h5ad output artifact
- `sc_tools/cli/benchmark.py` -- example command with multiple artifacts

### Clustering functions (PRV-05 targets)
- `sc_tools/bm/integration.py` -- `compute_integration_metrics()` calls `sc.tl.leiden`
- `sc_tools/pp/` -- preprocessing functions that may call clustering

### Pipeline and registry
- `sc_tools/pipeline.py` -- `STANDARD_PHASES`, phase DAG for lineage context

### Requirements
- `.planning/REQUIREMENTS.md` -- PRV-01 through PRV-05

### Prior context
- `.planning/phases/02-cli-foundation/02-CONTEXT.md` -- CLIResult design, error taxonomy, Provenance D-10
- `.planning/phases/03-core-commands/03-CONTEXT.md` -- Command structure, artifacts field usage

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `Provenance` Pydantic model in `models/result.py` -- expand with full fields (currently minimal)
- `cli_handler` decorator in `cli/__init__.py` -- wraps every command, has access to CLIResult after execution
- `CLIResult.artifacts` list -- already populated by commands with output file paths
- `_get_version()` in `models/result.py` -- version string already available

### Established Patterns
- Lazy imports inside command functions (CLI-08) -- new provenance commands must follow this
- `@cli_handler` decorator on every command function -- sidecar writing hooks into this
- `add_typer()` registration in `__init__.py` -- new `provenance` command group registers here
- `register_discovery(app)` pattern from Phase 4 -- follow for provenance command registration

### Integration Points
- `cli_handler` after `_emit()` -- insert sidecar writing logic
- `adata.uns['sct_provenance']` -- write provenance dict when output is h5ad
- `compute_integration_metrics(resolution, random_state)` -- thread random_state through
- New `cli/provenance.py` module for `show` and `trace` commands
- `models/result.py` `Provenance` model expansion (backwards compatible -- new fields optional for existing commands)

</code_context>

<specifics>
## Specific Ideas

- Hybrid provenance storage solves file-move portability: sidecar is canonical, adata.uns is portable backup. When files move between scratch/archive storage, SHA256 matching recovers the lineage chain.
- Relative paths (to project root) prevent absolute path breakage when projects move between machines or HPC scratch directories.
- All-flags-as-passed param serialization means agents can reconstruct the exact `sct` command from provenance alone — important for reproducibility.

</specifics>

<deferred>
## Deferred Ideas

- W3C PROV-JSON export (INF-01 in v2 requirements)
- RO-Crate packaging for publication (INF-02 in v2 requirements)
- Database-backed provenance migration (INF-03 in v2 requirements)
- DOT/graph visualization of lineage DAG -- start with flat list, add graph export if needed

</deferred>

---

*Phase: 05-provenance*
*Context gathered: 2026-03-23*
