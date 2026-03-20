# Codebase Concerns

**Analysis Date:** 2026-03-20

## Database Migration Debt

**Migrations 0012-0018 not applied to production:**
- Issue: Seven new migrations written (qc_status, four-layer schema, legacy table drops, samples, provenance, project_status_check) but not yet committed or applied
- Files: `sc_tools/migrations/versions/0012_add_qc_status.py`, `0013_four_layer_schema.py`, `0014_drop_legacy_tables.py`, `0015_drop_project_platform_domain.py`, `0016_project_status_check.py`, `0017_add_samples.py`, `0018_drop_provenance_inventory_id.py`
- Impact: Schema validation and registry operations may be out of sync between code and deployed database; four-layer schema redesign (separating data_sources, inventory_items, datasets) cannot be used until migrations run; patient/sample tracking incomplete
- Fix approach: Create migration control script, dry-run on test database, apply to production with backup, verify schema with `alembic current`

**Four-layer schema migration complexity:**
- Issue: Migration 0013 creates 7 new tables and migrates data from old schema; contains cross-join data migration logic that may fail on large registries
- Files: `sc_tools/migrations/versions/0013_four_layer_schema.py` (lines 39+)
- Impact: Long migration window; if data counts are high (>100k rows), migration could timeout; rollback would be slow
- Fix approach: Pre-migrate with sampling, test on production-like data volume, consider batch migration with checkpoints

**Old tables coexist with new:**
- Issue: Migration 0013 is additive-only; old `data_processing_phase`, `data_inventory` tables remain alongside new `inventory_items`, `datasets`, `data_sources` tables. Dual-write pattern used in `register_dataset()` (line 827 in registry.py)
- Files: `sc_tools/registry.py` (lines 827-850), `sc_tools/biodata.py` (register_biodata function)
- Impact: Dead code path in registry operations; ORM models may be inconsistent; queries could hit wrong table; cleanup via migrations 0014+ required
- Fix approach: Execute drop migrations (0014+) to remove legacy tables; add unit test to verify single-write after cleanup

---

## Tech Debt: Legacy Phase Nomenclature

**Dual phase-name support adds maintenance burden:**
- Issue: Codebase accepts both old p-codes (p1, p2, p3, p35, p4) and new semantic slugs (qc_filter, metadata_attach, preprocess, scoring, celltype_manual). Deprecation warning emitted but both paths functional
- Files: `sc_tools/validate.py` (lines 1-50, _LEGACY_PHASE_MAP), `scripts/validate_checkpoint.py`, `projects/create_project.sh` (Snakefile template)
- Impact: Users may use old codes in new scripts; maintenance burden of supporting both; confusion in logs and phase tracking
- Fix approach: Set firm deprecation deadline (e.g., 2026-06-01); emit DeprecationWarning systematically; migrate all projects to slug names; remove old code path in v2.0

---

## Missing CosMx Loader

**CosMx support incomplete:**
- Issue: CosMx modality listed in CLI help (`scripts/run_qc_report.py` line 606) but loader raises NotImplementedError
- Files: `sc_tools/ingest/loaders.py` (line 190, comment mentions cellid format but no implementation), `sc_tools/tests/test_ingest_cosmx.py` (test class `TestCosmxNotImplemented`)
- Impact: Projects cannot ingest CosMx data; cosmx_6k and cosmx_full_library projects deleted; only cosmx_1k remains; users see false promise in CLI
- Fix approach: Either implement loader (requires NanoString data parser) or remove from CLI choices; update Plan.md to reflect deprioritization

---

## Large/Complex Modules at Risk

**QC module file size and coupling:**
- Issue: `sc_tools/qc/plots.py` (1713 lines), `sc_tools/qc/report.py` (1528 lines), `sc_tools/tests/test_qc.py` (2054 lines) are among largest files in codebase. Tightly coupled to HTML templates and report generation
- Files: `sc_tools/qc/plots.py`, `sc_tools/qc/report.py`, `sc_tools/qc/report_utils.py`
- Impact: Hard to modify without affecting multiple functions; test coverage scattered across 2000+ lines of test code; regression risk high when changing QC metrics or plot logic
- Fix approach: Refactor plots into sub-modules by report type (pre_filter, post_filter, integration, celltyping); extract plot building to factory functions; add integration tests per report type

**Registry module state management:**
- Issue: `sc_tools/registry.py` (1811 lines) manages SQLAlchemy sessions via context manager `_session()` (line 413) but session lifecycle unclear in error paths
- Files: `sc_tools/registry.py` (lines 413-1000+)
- Impact: Session leaks possible on exception; transaction semantics unclear for nested calls; potential deadlocks on concurrent access
- Fix approach: Add explicit session cleanup tests; document transaction model; consider connection pooling review

---

## Soft Dependencies Not Validated at Install

**Optional packages missing validation:**
- Issue: Soft dependencies (scVI, Harmony, CytoVI, Plotly, pyyaml, etc.) imported at function call time with try/except but no install-time validation. Users may install sc-tools but lack required extras
- Files: `sc_tools/pp/integrate.py` (lines 90-96), `sc_tools/qc/report_utils.py` (lines 167-169, 326, 403, 424), `sc_tools/data/imc/benchmark/config.py` (lines 94-95)
- Impact: Errors surface only when calling functions, not at import time; poor user experience if extras not installed; no guidance in error messages about which extras to install
- Fix approach: Add runtime check in `__init__.py` to warn about missing extras; improve error messages with install commands; consider optional test suite that skips if extras missing

**Plotly tests skip gracefully:**
- Issue: QC report tests skip when plotly missing (recent fix from commit 41326d5). Good coverage for plotly, but other optional deps not similarly guarded
- Files: `sc_tools/qc/report_utils.py` (line 167, ImportError handler), test files with skipif decorators
- Impact: Integration tests may pass locally but fail in CI if optional deps missing; inconsistent test behavior across environments
- Fix approach: Add CI matrix for optional/no-optional installs; document which tests require which extras

---

## Unimplemented Features Stubbed

**Multiple backend stubs in celltype module:**
- Issue: `sc_tools/tl/celltype/` has stubs for scArches, Geneformer, scGPT, SingleR that raise NotImplementedError
- Files: `sc_tools/tl/celltype/_scarches.py`, `_geneformer.py`, `_scgpt.py`, `_singler.py` (all line 20, `raise NotImplementedError(...)`)
- Impact: Users may discover unimplemented methods only at runtime; test suite still validates stubs but does not implement functionality; dead code in repo
- Fix approach: Implement or remove stubs; if deferred, add clear documentation of timeline and blockers; remove from public API until ready

**Geary C spatial autocorrelation not supported:**
- Issue: `sc_tools/gr/_autocorr.py` (line 50, 64) raises NotImplementedError for Geary C (only Moran I supported)
- Files: `sc_tools/gr/_autocorr.py`, test validates rejection at line 323 `test_gr.py`
- Impact: Users cannot compute Geary C; API suggests support but fails at runtime; design choice not documented
- Fix approach: Document squidpy constraint; add Geary C support if demand exists; or remove from API entirely

**resolVI integration deprioritized:**
- Issue: `run_resolvi()` listed in Plan.md as future work; not implemented
- Files: `sc_tools/pp/integrate.py` (line 28, in __all__)
- Impact: Function stub may mislead users; no timeline for implementation
- Fix approach: Remove from exports if not planned; or add timeline to docstring

---

## Test Coverage Gaps

**Real-data tests may be stale:**
- Issue: Real-data test files exist (`test_qc_real_data.py`, `test_pp_real_data.py`, `test_tl_real_data.py`, `test_gr_real_data.py`) but depend on project data that may have moved or been deleted
- Files: `sc_tools/tests/test_*_real_data.py`
- Impact: Real-data tests may fail silently if data paths change; no CI validation of test data freshness; coverage for edge cases (empty ROI, sparse cell types) incomplete
- Fix approach: Standardize test data paths using fixtures; validate data freshness in CI; add skip-if-no-data logic

**Migration test gaps:**
- Issue: Migration tests skipped without alembic (test_migrations.py line 118, `pass` statement)
- Files: `sc_tools/tests/test_migrations.py`
- Impact: Migration logic never validated in CI; schema changes may break silently; rollback paths untested
- Fix approach: Always include alembic in test dependencies; add forward/backward migration tests; validate schema integrity after each migration

**Edge case coverage incomplete:**
- Issue: Plan.md (line 121) lists "Edge-case coverage: empty ROI, sparse cell types, spatial grid subsets, multi-library_id" as TODO
- Files: Multiple qc, pp, gr test files
- Impact: Error handling for edge cases may be missing; crashes possible on empty data or unusual sample configurations
- Fix approach: Add edge-case test factory; parametrize tests with empty/sparse data; verify all modalities handle edge cases

---

## Integration Points Fragile

**HTML template and report code tight coupling:**
- Issue: Report generation code in `sc_tools/qc/report.py` is tightly coupled to Jinja2 templates in `sc_tools/assets/*.html`. Changes to template structure require code changes and vice versa. Tabbed navigation (Option B) uses regex to parse prior HTML reports
- Files: `sc_tools/qc/report.py` (lines 1-100), `sc_tools/qc/report_utils.py` (lines 4-41, _wrap_with_tabs, _extract_body_content, _extract_head_css), `sc_tools/assets/pre_filter_qc_template.html`, `post_filter_qc_template.html`
- Impact: Template changes may silently break reports; regex-based HTML parsing fragile to whitespace changes; no validation that generated HTML is well-formed
- Fix approach: Use proper HTML parser (BeautifulSoup) instead of regex; add integration tests that generate HTML and validate structure; separate template context builders from rendering

**Benchmark integration config complexity:**
- Issue: `sc_tools/bm/integration.py` (426 lines) has conditional branching for integration_method selection and parameter passing. TODO comment at line 426 ("When unintegrated baseline wins, consider emitting a warning") suggests unresolved design
- Files: `sc_tools/bm/integration.py` (lines 1-426)
- Impact: Code path coverage unclear; parameter validation missing; may select wrong method under edge conditions
- Fix approach: Add comprehensive tests for each method; add validation of method-specific params; document method selection logic

---

## Performance Bottlenecks

**K-NN subsampling for large datasets:**
- Issue: `sc_tools/qc/plots.py` (lines 1403-1431, k-NN contextual pass) subsamples to 50k points for fitting NearestNeighbors, which may impact accuracy for larger samples
- Files: `sc_tools/qc/plots.py` (line 1403, FIT_LIMIT = 50_000)
- Impact: High-res spatial data (>50k spots) loses neighborhood context precision; false negatives in viable region assessment possible
- Fix approach: Benchmark accuracy vs performance; consider hierarchical clustering for large datasets; or increase FIT_LIMIT with memory profiling

**TPM metrics computation inefficiency:**
- Issue: `compute_sample_metrics()` applies viable_threshold and TPM calculations in-memory per sample, iterating through all samples. No vectorization
- Files: `sc_tools/qc/sample_qc.py` (lines 200-300 approx, compute_sample_metrics signature)
- Impact: Slow for projects with 100+ samples; no caching between calls
- Fix approach: Vectorize computations using pandas/numpy operations; cache intermediate results; profile on large datasets

---

## Security/Validation Gaps

**URI validation weak:**
- Issue: `sc_tools/registry.py` (_validate_uri, line 97-100) checks if URI starts with known schemes but does not validate path existence, permissions, or sanitize for injection
- Files: `sc_tools/registry.py` (lines 97-105)
- Impact: Invalid URIs accepted silently; potential for path traversal if used in shell commands; no feedback if data source becomes unreachable
- Fix approach: Add URI scheme-specific validation; test path accessibility at registration time; warn on stale URIs

**Clinical data injection risk:**
- Issue: `scripts/seed_clinical_data.py` parses XLSX files and inserts into database. No validation of data types or ranges; potential for SQL injection if ORM not used correctly
- Files: `scripts/seed_clinical_data.py` (lines 1-600)
- Impact: Malformed data could corrupt registry; user input not sanitized
- Fix approach: Add schema validation layer; type-check all fields; add unit tests for invalid XLSX inputs

**Environment variable leakage:**
- Issue: Registry uses SC_TOOLS_REGISTRY_URL env var for DB connection. If DB URL contains credentials, could be logged or leaked in error messages
- Files: `sc_tools/registry.py` (lines 69-76, _get_db_url)
- Impact: Accidental credential exposure in logs or error output
- Fix approach: Mask credentials in debug output; use separate auth env vars; validate URL does not contain passwords

---

## Documentation and Clarity Gaps

**Phase nomenclature confusion:**
- Issue: Code uses both semantic slugs (qc_filter, metadata_attach) and old p-codes (p1, p2, p3, p35, p4). Documentation in Architecture.md and Plan.md inconsistent
- Files: `docs/Architecture.md`, `docs/Plan.md`, `sc_tools/validate.py` (line 35-47)
- Impact: New users confused about phase naming; commits use both notations; migration to new scheme unclear
- Fix approach: Unified documentation with mapping table; enforce slug usage in new code; deprecate p-codes by 2026-06-01

**Viable region assessment under-documented:**
- Issue: Viable region assessment (added 2026-03-18) uses viable_threshold parameters (min_counts, min_genes) but semantics differ from QC filter thresholds. Not clear to users when to adjust
- Files: `docs/Journal.md` (lines 27-37), `sc_tools/qc/sample_qc.py` (viable_threshold param)
- Impact: Users may misconfigure thresholds; report interpretations unclear; no guidance on modality-specific tuning
- Fix approach: Add section to Architecture.md explaining viable region logic; add docstring examples for typical thresholds by modality

**Biodata hierarchy complexity:**
- Issue: BioData v1.1 introduces modality tier, PlatformSpec defaults, clinical data dual-write (Journal 2026-03-11). 3-tier classification (BioDataType -> BioDataModality -> BioDataPlatform) documented only in Journal.md, not in main Architecture.md
- Files: `docs/Journal.md` (lines 39-54), `sc_tools/biodata.py`
- Impact: Developers unfamiliar with biodata structure; unclear how to extend taxonomy; dual-write logic not obvious
- Fix approach: Add comprehensive biodata section to Architecture.md; create biodata_api.md reference guide; docstring examples for register_biodata()

---

## Inconsistent Error Handling

**ImportError messages vary:**
- Issue: Soft dependency ImportError messages inconsistent in formatting and detail. Some provide install command, others do not
- Files: `sc_tools/pp/integrate.py` (line 93-96), `sc_tools/qc/report_utils.py` (line 168-169), `sc_tools/storage.py` (line 70, 80)
- Impact: Users get different quality of error messages; some may not know how to fix import errors
- Fix approach: Create unified error message template; apply consistently across all soft dependencies; include install command in every ImportError

**Session context manager not documented:**
- Issue: Registry uses `with self._session()` pattern (line 413+) but cleanup behavior on exception unclear
- Files: `sc_tools/registry.py` (lines 413-450)
- Impact: Potential for resource leaks; error recovery path unknown
- Fix approach: Document context manager contract; add explicit finally block; write tests for exception handling

---

## Scaling Limits

**Registry SQLite default:**
- Issue: Registry defaults to SQLite at ~/.sc_tools/registry.db (line 61, _DEFAULT_DB_PATH). SQLite has write locks; not suitable for concurrent access
- Files: `sc_tools/registry.py` (lines 61, 75-76), `sc_tools/mcp/registry_server.py` (MCP server may have concurrent requests)
- Impact: High contention environments (shared HPC, multi-agent) will see lock timeouts; data corruption possible under write load; Postgres option exists but not default
- Fix approach: Document Postgres requirement for production; add connection pooling if using SQLite; add locking tests for concurrent scenarios

**Array storage in JSON metadata:**
- Issue: `data_sources`, `inventory_items`, `datasets` tables use JSONB/JSON metadata columns (migrations 0013). Large arrays in metadata not handled efficiently
- Files: `sc_tools/migrations/versions/0013_four_layer_schema.py` (lines 60, 99, etc., metadata columns)
- Impact: JSON queries slow; partial indexes needed for common metadata fields; schema not normalized for array access
- Fix approach: Normalize frequently-queried metadata into separate columns; add indexes; profile query performance on real data

---

## Known Bugs Fixed Recently

**K-NN subsampling index mismatch (fixed 2026-03-18b):**
- Issue: Bug 1 in Journal (2026-03-18b): Index mismatch in k-NN subsampling path; code used `_indices` directly without mapping back through subsampled indices
- Files: `sc_tools/qc/sample_qc.py`, `sc_tools/qc/plots.py` (line 1426, now fixed: `tc[sub_idx[_indices]]`)
- Impact: Viable region assessment gave wrong neighbor counts for large samples; false pass/fail classifications possible
- Status: Fixed in commit; requires verification on real 60K+ spot samples

**Template column header confusion (fixed 2026-03-18b):**
- Issue: Bug 3: Pre/post filter templates showed "Stromal Cells" instead of "Stromal Spots" for Visium; used linear counts in plot (2,1) instead of log
- Files: `sc_tools/assets/pre_filter_qc_template.html`, `post_filter_qc_template.html` (now use `terms.observations_lower`)
- Status: Fixed; template headers now dynamic based on modality

**QC 2x2 grid panel layout (fixed 2026-03-18b):**
- Issue: Bug 5: Panel (1,1) showed raw counts instead of log10; panel (2,2) logic unclear
- Files: `sc_tools/qc/plots.py` (now log10(counts+1) in panel 1,1)
- Status: Fixed; panel logic now clearer with conditional MT/genes display

---

## Remaining Incomplete Features

**Post-normalization QC report:**
- Issue: Plan.md (line 76) lists "Post-normalization QC report" as TODO under preprocess phase
- Files: None yet
- Impact: Users cannot validate normalization quality; normalization parameters uncheckable
- Fix approach: Design report with before/after normalization metrics; add to preprocess phase

**Demographics figure generation:**
- Issue: Plan.md (lines 80) lists "sc_tools helpers: piechart, histogram, violin, bar, heatmap" and "Figure 1 for cohort description" as TODO
- Files: None yet
- Impact: Cannot generate cohort demographic figures; Figure 1 must be done manually
- Fix approach: Implement plotting helpers; add demographics report template

**Project-local Snakemake workflows for IMC:**
- Issue: Plan.md (line 37) lists "IMC end-to-end Snakemake workflow for lymph_dlbcl and ggo-imc" as TODO; ingest_load Snakemake rule missing (line 38)
- Files: Projects `projects/imc/lymph_dlbcl/`, `projects/imc/ggo_imc/` have no Snakefiles
- Impact: IMC projects cannot be run via pipeline; manual script execution only
- Fix approach: Implement IMC Snakemake workflow; add ingest_load sentinel

---

## Workarounds and Provisional Decisions

**Phase B dual-write deprecation path:**
- Issue: `register_dataset()` writes to both old `data_processing_phase` and new schema. Phase A (warning) -> Phase B (dual-write) -> Phase C (removal) planned (Journal, line 47)
- Files: `sc_tools/biodata.py` (register_biodata function)
- Impact: Dead code path; confusion about which tables to query; cleanup deferred
- Fix approach: Execute Phase C removal when migrations 0014+ are applied; add deprecation deadline

**Conditional template selection in reports:**
- Issue: `generate_post_celltyping_report()` was using wrong template (_POST_INTEGRATION_TEMPLATE instead of _POST_CELLTYPING_TEMPLATE, Journal line 76). Fixed by building ranking_rows explicitly
- Files: `sc_tools/qc/report.py` (line 59, now correctly references _POST_CELLTYPING_TEMPLATE)
- Status: Fixed; but indicates template selection logic fragile

---

## Summary of Priority Concerns

| Concern | Severity | Blocker | Timeline |
|---------|----------|---------|----------|
| Migrations 0012-0018 not applied | HIGH | Yes - schema mismatch | Apply immediately |
| CosMx loader unimplemented | MEDIUM | No - deprioritized | Remove or implement Q2 2026 |
| QC module size + coupling | MEDIUM | No - works but fragile | Refactor Q2 2026 |
| Registry SQLite concurrency | MEDIUM | Yes - under load | Switch to Postgres or add locking |
| Soft dependency validation | MEDIUM | No - UX issue | Add validation checks Q1 2026 |
| Plotly/optional test coverage | MEDIUM | No | Add CI matrix with/without optionals |
| Real-data test staleness | MEDIUM | No | Standardize test data fixtures Q2 2026 |
| Phase nomenclature debt | LOW | No | Enforce slugs, deprecate p-codes by 2026-06-01 |
| Documentation gaps (biodata, viable regions) | LOW | No | Update Architecture.md Q1 2026 |

