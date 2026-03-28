# Phase 9: Sample Concatenation & Maintenance - Research

**Researched:** 2026-03-27
**Domain:** AnnData concatenation CLI, pipeline DAG registration, Plotly dependency management
**Confidence:** HIGH

## Summary

Phase 9 has two independent workstreams: (1) a new `sct concat` CLI command that wraps AnnData concatenation with spatial data preservation and provenance tracking, and (2) two maintenance fixes to the Plotly dependency chain. Both are well-scoped changes that follow established codebase patterns.

The concat command is primarily a CLI wrapper around `anndata.concat()` with the critical parameter `uns_merge="unique"` to preserve per-sample spatial metadata. The existing `concat_samples()` function in `ingest/loaders.py` uses `uns_merge="same"`, which empirically **drops all spatial keys** when samples have different library IDs (verified via test). The new CLI command must either fix this function or bypass it to use `ad.concat` directly with the correct merge strategy.

The Plotly maintenance requires two changes: promoting `plotly>=5.18` from the `[pipeline]` optional extra to base `dependencies` in `pyproject.toml`, and updating the CDN pin from `plotly-2.27.0.min.js` to `plotly-3.4.0.min.js`. The REQUIREMENTS.md specifies replacing with `plotly-latest.min.js`, but this is incorrect -- `plotly-latest.min.js` is frozen at v1.58.5 on the CDN and will never be updated. The installed plotly.py 6.6.0 bundles plotly.js 3.4.0, so the CDN must match. The project's frontend skill (`k-dense-frontend/SKILL.md`) explicitly warns: "NEVER use `plotly-latest` -- the CDN deprecated it and it may break silently."

**Primary recommendation:** Implement as 3 plans: (1) Plotly dependency promotion + CDN update, (2) `sct concat` CLI command with `uns_merge="unique"` + provenance sidecar, (3) pipeline DAG registration of `concat` phase.

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| MAINT-01 | `plotly` promoted from `[pipeline]` optional to base `dependencies` | Move `"plotly>=5.18"` from `[pipeline]` list to `[project.dependencies]` in pyproject.toml. Also remove from `[benchmark]` list to avoid duplication. Verified: plotly 6.6.0 installed, works with Python 3.10+. |
| MAINT-02 | Stale CDN pin replaced with current version | Replace `plotly-2.27.0.min.js` with `plotly-3.4.0.min.js` in 3 locations: `report_utils.py` line 517, `base_report_template.html` line 8, and the frontend skill SKILL.md. **Do NOT use `plotly-latest.min.js`** -- it is frozen at v1.58.5. Update test assertion in `test_qc.py` line 1311. |
| CONCAT-01 | `sct concat` with list of h5ad paths produces merged h5ad | New CLI command in `sc_tools/cli/concat.py` following `register_concat(app)` pattern (like `assemble.py`). Accepts `--input` list + `--output` path. Uses `anndata.concat()` with `join="outer"`, `uns_merge="unique"`, `index_unique="-"`. |
| CONCAT-02 | Preserves `uns["spatial"]` via `uns_merge="unique"` | Empirically verified: `uns_merge="unique"` preserves all per-sample spatial keys. `uns_merge="same"` (current `concat_samples` default) drops them all. Post-concat verification: assert input library IDs == output `uns['spatial']` keys. |
| CONCAT-03 | Writes `.provenance.json` sidecar with SHA256 checksums | Use existing `_input_files` convention in CLIResult.data -- the `cli_handler` decorator auto-calls `build_provenance_record()` + `write_sidecar()` which computes SHA256 via `provenance/checksum.py`. No new code needed for sidecar writing. |
| CONCAT-04 | Registered as pipeline phase `concat` between `ingest_load` and `qc_filter` | Add `"concat"` PhaseSpec to `STANDARD_PHASES` in `pipeline.py` with `depends_on=[_dp("ingest_load")]`, `optional=True`. Update `qc_filter` depends_on to include concat OR keep as-is since optional phases can be skipped. |
</phase_requirements>

## Standard Stack

### Core

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| anndata | 0.11.4 | `ad.concat()` for merging h5ad files | Already installed; native merge with `uns_merge` parameter |
| plotly | 6.6.0 (py) / 3.4.0 (js) | Interactive charts in HTML reports | Already installed; bundles plotly.js 3.4.0 |
| typer | 0.24.1 | CLI framework | Already used for all `sct` commands |
| pydantic | 2.10+ | CLIResult envelope | Already used for structured CLI output |
| h5py | (bundled) | Pre-flight validation without full load | Read obs/var index from h5ad without loading X matrix |

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| hashlib | stdlib | SHA256 checksums for provenance | Already used in `provenance/checksum.py` |
| scanpy | 1.9+ | `sc.pp.calculate_qc_metrics` post-concat | Optional QC recalculation after merge |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `uns_merge="unique"` | `uns_merge="first"` | Both preserve spatial keys empirically, but "unique" is semantically correct -- each sample's spatial dict is unique |
| New `cli/concat.py` | Extending `cli/qc.py` | Concat is logically an ingestion operation, not QC -- separate module is cleaner |
| Fixing `concat_samples()` | Bypass it entirely | Fixing the existing function is better for backward compat, but CLI should validate result regardless |

## Architecture Patterns

### Recommended Project Structure (changes only)

```
sc_tools/
  cli/
    concat.py          # NEW - sct concat command (register_concat pattern)
    __init__.py         # MODIFIED - register concat_app
  pipeline.py           # MODIFIED - add "concat" PhaseSpec to STANDARD_PHASES
  ingest/
    loaders.py          # MODIFIED - fix concat_samples uns_merge="same" -> "unique"
  qc/
    report_utils.py     # MODIFIED - CDN pin update
  assets/
    base_report_template.html  # MODIFIED - CDN pin update
  tests/
    test_concat.py      # NEW - concat CLI and logic tests
    test_qc.py          # MODIFIED - update CDN assertion
    test_pipeline.py    # MODIFIED - test concat phase in DAG
```

### Pattern 1: CLI Command Registration (established pattern from assemble.py)

**What:** New CLI commands use the `register_X(app)` pattern with a nested Typer app.
**When to use:** Any new top-level `sct` subcommand.
**Example:**
```python
# Source: sc_tools/cli/assemble.py (existing pattern)
def register_concat(app: typer.Typer) -> None:
    """Register the concat command on the given Typer app."""
    from sc_tools.cli import _check_deps, cli_handler
    from sc_tools.models.result import CLIResult, Provenance, Status

    @app.command("concat")
    @cli_handler
    def concat(
        input: list[str] = typer.Option(..., "--input", "-i", help="Input h5ad files"),
        output: str = typer.Option(..., "--output", "-o", help="Output merged h5ad"),
        batch_key: str = typer.Option("sample", "--batch-key", help="Obs column for sample identity"),
        dry_run: bool = typer.Option(False, "--dry-run", help="Validate without merging"),
        force: bool = typer.Option(False, "--force", help="Bypass memory guard"),
    ) -> CLIResult:
        ...
```

### Pattern 2: Provenance Sidecar via _input_files Convention

**What:** CLI commands include `_input_files` in `CLIResult.data` and the `cli_handler` decorator auto-generates the sidecar.
**When to use:** Any command that produces output artifacts.
**Example:**
```python
# Source: sc_tools/cli/__init__.py lines 283-296 (existing pattern)
return CLIResult(
    status=Status.success,
    command="concat",
    data={
        "n_samples": len(adatas),
        "total_cells": merged.n_obs,
        "spatial_keys_preserved": list(merged.uns.get("spatial", {}).keys()),
        "_input_files": input_paths,  # triggers automatic sidecar
    },
    artifacts=[output],
    provenance=Provenance(command="concat"),
    message=f"Concatenated {len(adatas)} samples: {merged.n_obs} cells",
)
```

### Pattern 3: Pipeline Phase Registration

**What:** Add a PhaseSpec to `STANDARD_PHASES` dict for DAG visibility.
**When to use:** New recognized pipeline phases.
**Example:**
```python
# Insert between "ingest_load" and "qc_filter" in pipeline.py STANDARD_PHASES
"concat": PhaseSpec(
    label="Sample Concatenation",
    depends_on=[_dp("ingest_load")],
    branch="ingestion",
    checkpoint="results/adata.concatenated.h5ad",
    phase_group=_DP,
    required_obs=["sample"],
    required_obsm=["spatial"],
    x_format="raw counts, concatenated",
    optional=True,
),
```

### Anti-Patterns to Avoid

- **Using `uns_merge="same"` for concat:** Drops all spatial data when samples have different library IDs. Always use `uns_merge="unique"`.
- **Loading full AnnData for pre-flight validation:** Use h5py to read only obs/var indices. Do NOT `ad.read_h5ad()` just to check compatibility.
- **Using `plotly-latest.min.js` CDN:** Frozen at v1.58.5 (July 2021). Always pin an explicit version.
- **Skipping `.copy()` on AnnData subset before spatial plots:** View slices may not carry `uns["spatial"]` correctly.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| SHA256 checksums | Custom hashing | `sc_tools.provenance.checksum.sha256_file()` | Already implemented with streaming reads |
| Provenance sidecar | Manual JSON writing | `cli_handler` decorator + `_input_files` convention | Automatic, atomic writes, h5ad embedding |
| CLI error handling | Try/except in command | `@cli_handler` decorator | Maps exceptions to exit codes, emits CLIResult |
| Phase DAG queries | Manual dependency resolution | `pipeline.get_available_next()` | Handles optional, iterative, and custom phases |

**Key insight:** The entire provenance pipeline is automatic when you follow the `_input_files` convention. The `cli_handler` decorator intercepts the CLIResult, calls `build_provenance_record()`, and writes the sidecar atomically.

## Common Pitfalls

### Pitfall 1: uns["spatial"] Loss Through Concat
**What goes wrong:** After `ad.concat()`, spatial keys are missing and downstream `sc.pl.spatial()` calls fail with KeyError.
**Why it happens:** `uns_merge="same"` (current default in `concat_samples`) only keeps uns entries that are identical across all inputs. Since each sample has a unique spatial dict, all are dropped.
**How to avoid:** Use `uns_merge="unique"` and post-concat verification: `assert set(expected_lib_ids) == set(merged.uns.get('spatial', {}).keys())`
**Warning signs:** `merged.uns` is empty or missing `spatial` key after concat.

### Pitfall 2: plotly-latest.min.js is v1.58.5
**What goes wrong:** CDN loads plotly.js v1.x instead of v3.x, and charts generated by plotly.py 6.x fail to render or display incorrectly.
**Why it happens:** Plotly stopped updating `plotly-latest.min.js` at v1.58.5 when v2.0 shipped (breaking API change). The filename is misleading.
**How to avoid:** Always pin an explicit version. Currently `plotly-3.4.0.min.js` matches the bundled plotly.js in plotly.py 6.6.0.
**Warning signs:** Charts silently degrade or JS console errors about unknown trace types.

### Pitfall 3: Duplicate obs Index After Concat
**What goes wrong:** Multiple samples share cell barcodes (e.g., "AACG..."), causing duplicate index entries that break downstream operations.
**Why it happens:** 10x Genomics barcodes are shared across samples.
**How to avoid:** Use `index_unique="-"` in `ad.concat()` to append sample key to index.
**Warning signs:** `merged.obs_names_are_unique` returns False.

### Pitfall 4: Modality Mismatch in Concat Inputs
**What goes wrong:** User passes h5ad files with different var_names (e.g., one Visium and one IMC sample), producing a sparse matrix with mostly zeros.
**Why it happens:** `join="outer"` fills missing genes with zeros/NaN.
**How to avoid:** Pre-flight check: read var_names from each input via h5py and warn if overlap < 80%.
**Warning signs:** Output h5ad has far more vars than any single input.

### Pitfall 5: Test Assertion Hardcodes CDN Version
**What goes wrong:** Tests pass before CDN update but fail after.
**Why it happens:** `test_qc.py` line 1311 asserts `"plotly-2.27.0.min.js" in content`.
**How to avoid:** Update the test assertion alongside the CDN pin change.
**Warning signs:** CI failure in test_qc after CDN update.

## Code Examples

### Example 1: Correct AnnData Concat with Spatial Preservation

```python
# Verified empirically on anndata 0.11.4
import anndata as ad

merged = ad.concat(
    adatas,
    join="outer",
    uns_merge="unique",   # preserves per-sample spatial dicts
    index_unique="-",     # deduplicates cell barcodes
    label="sample",       # adds sample column to obs
    keys=sample_names,    # sample identifiers
)

# Post-concat verification (REQUIRED)
expected = {name for name in sample_names}
actual = set(merged.uns.get("spatial", {}).keys())
missing = expected - actual
if missing:
    raise SCToolsDataError(
        f"Spatial data lost for samples: {missing}",
        suggestion="Check input h5ad files have uns['spatial'] populated",
    )
```

### Example 2: h5py Pre-Flight Validation (No Full Load)

```python
import h5py

def validate_h5ad_preflight(path: str) -> dict:
    """Read metadata from h5ad without loading X matrix."""
    with h5py.File(path, "r") as f:
        n_obs = f["obs"].attrs.get("_index", f["obs/_index"]).shape[0]
        var_names = f["var/_index"][()].astype(str).tolist()
        has_spatial = "spatial" in f.get("uns", {})
    return {"n_obs": n_obs, "var_names": var_names, "has_spatial": has_spatial}
```

### Example 3: Pipeline PhaseSpec for Optional Concat

```python
# In pipeline.py STANDARD_PHASES, between ingest_load and qc_filter
"concat": PhaseSpec(
    label="Sample Concatenation",
    depends_on=[_dp("ingest_load")],
    branch="ingestion",
    checkpoint="results/adata.concatenated.h5ad",
    phase_group=_DP,
    required_obs=["sample"],
    required_obsm=["spatial"],
    x_format="raw counts, concatenated",
    optional=True,  # projects can skip if pre-concatenated
),
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `plotly-2.27.0.min.js` CDN | `plotly-3.4.0.min.js` CDN | plotly.py 6.0 (Jan 2025) | Charts use plotly.js v3 API; CDN must match |
| `plotly-latest.min.js` | Explicit version pin | plotly.js v2.0 (2021) | `-latest` frozen at v1.58.5, never updated |
| `uns_merge="same"` in concat | `uns_merge="unique"` | Always was the correct option | "same" drops spatial data; "unique" preserves it |
| plotly in `[pipeline]` extra | plotly in base deps | v2.0 (this phase) | Reports available without `pip install .[pipeline]` |

**Deprecated/outdated:**
- `plotly-latest.min.js`: Frozen at v1.58.5 since plotly.js v2.0 release. Never use.
- `concat_samples()` with `uns_merge="same"`: Loses spatial data. Must be fixed.

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | pytest 7.0+ |
| Config file | `pyproject.toml` [tool.pytest.ini_options] |
| Quick run command | `pytest sc_tools/tests/test_concat.py -x -q` |
| Full suite command | `pytest sc_tools/tests/ -q` |

### Phase Requirements -> Test Map

| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| MAINT-01 | plotly importable without `[pipeline]` | unit | `python -c "import plotly"` | N/A (pyproject.toml change) |
| MAINT-02 | CDN pin updated to 3.4.0 | unit | `pytest sc_tools/tests/test_qc.py::TestBaseTemplate -x` | Exists (needs assertion update) |
| CONCAT-01 | `sct concat` merges h5ad files | unit | `pytest sc_tools/tests/test_concat.py::TestConcatCommand -x` | Wave 0 |
| CONCAT-02 | `uns["spatial"]` preserved | unit | `pytest sc_tools/tests/test_concat.py::TestSpatialPreservation -x` | Wave 0 |
| CONCAT-03 | provenance sidecar with SHA256 | unit | `pytest sc_tools/tests/test_concat.py::TestConcatProvenance -x` | Wave 0 |
| CONCAT-04 | concat in pipeline DAG | unit | `pytest sc_tools/tests/test_pipeline.py -x -k concat` | Wave 0 |

### Sampling Rate
- **Per task commit:** `pytest sc_tools/tests/test_concat.py -x -q && make lint`
- **Per wave merge:** `pytest sc_tools/tests/ -q`
- **Phase gate:** Full suite green before `/gsd:verify-work`

### Wave 0 Gaps
- [ ] `sc_tools/tests/test_concat.py` -- covers CONCAT-01, CONCAT-02, CONCAT-03
- [ ] Pipeline test additions in `test_pipeline.py` -- covers CONCAT-04
- [ ] Update assertion in `test_qc.py` line 1311 -- covers MAINT-02

## Open Questions

1. **Should `concat_samples()` in `loaders.py` be fixed or left as-is?**
   - What we know: It currently uses `uns_merge="same"` which drops spatial data. It is called in the existing ingest pipeline.
   - What's unclear: Whether any existing pipeline code depends on the `uns_merge="same"` behavior.
   - Recommendation: Fix it to use `uns_merge="unique"` since "same" silently drops data. No downstream code should depend on spatial data being absent.

2. **Should `qc_filter` depend on `concat` in the DAG?**
   - What we know: `concat` is `optional=True`, so `qc_filter` cannot require it.
   - What's unclear: Whether `qc_filter.depends_on` should list both `ingest_load` and `concat` (with concat optional).
   - Recommendation: Keep `qc_filter.depends_on = [_dp("ingest_load")]` unchanged. The concat phase is optional and its output replaces the qc_filter input, but does not gate it. Projects that skip concat feed per-sample h5ads directly to qc_filter.

3. **MAINT-02 wording says `plotly-latest.min.js` -- should we follow it literally?**
   - What we know: `plotly-latest.min.js` is frozen at v1.58.5 and will break charts generated by plotly.py 6.x.
   - Recommendation: Override MAINT-02 wording. Use `plotly-3.4.0.min.js` (matches installed plotly.py 6.6.0). Document deviation clearly.

## Project Constraints (from CLAUDE.md)

- **Linting:** Ruff; never commit failing lint. Run `make lint` before every commit.
- **Imports:** scanpy as `sc`, anndata as `ad`, pandas as `pd`, numpy as `np`. Import sc_tools functions from `sc_tools.<module>`.
- **Checkpoint I/O:** Always use `sc_tools.io.write_checkpoint()` / `read_checkpoint()` -- never raw `adata.write()`.
- **Logging:** Use `sc_tools.utils.get_logger(__name__)`, not `print()`.
- **Tests:** Unit + integration; fail-proof with empty/sub/full fixtures.
- **CLI pattern:** `@cli_handler` decorator + `CLIResult` envelope for all commands.
- **Frontend reports:** Bootstrap 5.3 Flatly + Plotly (pinned version, NEVER `plotly-latest`). Base64 PNG for spatial plots.
- **TDD order:** Write failing test -> confirm RED -> implement -> confirm GREEN -> export -> lint -> commit.
- **No heavy imports at module level:** CLI-08 compliance (lazy import scanpy, h5py, etc.).

## Sources

### Primary (HIGH confidence)
- **Codebase inspection:** `pipeline.py`, `cli/__init__.py`, `cli/assemble.py`, `provenance/sidecar.py`, `ingest/loaders.py`, `qc/report_utils.py`, `models/result.py` -- all patterns verified by reading source
- **Empirical test:** `ad.concat()` with `uns_merge="unique"` vs `"same"` vs `"first"` -- verified on anndata 0.11.4
- **Empirical test:** `plotly.offline.get_plotlyjs_version()` returns `3.4.0` for plotly.py 6.6.0
- **Frontend skill:** `.claude/skills/k-dense-frontend/SKILL.md` -- "NEVER use `plotly-latest`"

### Secondary (MEDIUM confidence)
- [plotly.js CDN issue #7315](https://github.com/plotly/plotly.js/issues/7315) -- confirms `plotly-latest.min.js` frozen at v1.58.5
- [plotly.js releases](https://github.com/plotly/plotly.js/releases) -- current version is 3.4.0
- [plotly.py v6 migration](https://plotly.com/python/v6-migration/) -- plotly.py 6.x uses plotly.js 3.x

### Tertiary (LOW confidence)
- None -- all findings verified via codebase inspection or empirical testing.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- all libraries already installed and verified; no new dependencies
- Architecture: HIGH -- all patterns directly observed in existing codebase (assemble.py, cli_handler, provenance)
- Pitfalls: HIGH -- uns_merge behavior empirically verified; plotly CDN issue confirmed via GitHub issue

**Research date:** 2026-03-27
**Valid until:** 2026-04-27 (stable -- no fast-moving dependencies)
