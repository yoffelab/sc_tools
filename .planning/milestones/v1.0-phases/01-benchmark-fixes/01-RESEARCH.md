# Phase 1: Benchmark Fixes - Research

**Researched:** 2026-03-20
**Domain:** Integration benchmarking (scib-metrics, h5py, AnnData, NumPy)
**Confidence:** HIGH

## Summary

Phase 1 is a fix-and-harden phase for `sc_tools/bm/integration.py` and `sc_tools/pp/recipes.py`. All seven BM requirements and three TST requirements target existing code with well-understood bugs. The codebase already has 32 integration benchmark tests in `test_integration_benchmark.py` and a comprehensive test pattern to follow.

The primary technical challenges are: (1) reading obsm embeddings and obs columns from h5ad files via h5py without loading full AnnData (BM-01), (2) per-embedding NaN masking before metric computation (BM-02), and (3) fixing a truncation bias in `_stratified_subsample` (BM-04). The remaining requirements are straightforward additions (timing, provenance columns, conditional normalization skip).

**Primary recommendation:** Fix the seven bugs in `integration.py` and `recipes.py` with TDD (red-green-refactor), using synthetic AnnData fixtures. No new dependencies needed -- h5py is already installed (v3.15.1).

<user_constraints>

## User Constraints (from CONTEXT.md)

### Locked Decisions
- Add `embedding_files` parameter to `compare_integrations()`: dict mapping method name to h5ad file path
- Backwards compatible -- existing `adata` + `embeddings` dict API remains functional
- When `embedding_files` is provided, load only `obsm` embeddings + required `obs` columns via h5py (no full AnnData load)
- Extract `batch_key` plus a single optional `bio_key` (replaces the celltype-only assumption)
- `bio_key` defaults to `celltype_key` but can be any clinically relevant variable (e.g., `condition`, `disease_status`) for bio conservation metrics
- Cell type is often missing in early pipeline stages -- the API must work without it (batch-only metrics)
- Peak memory target: <2GB for 2.5M-cell dataset when using h5py path

### Claude's Discretion
- NaN handling strategy for resolVI embeddings (BM-02) -- per-embedding valid-cell masking is the natural approach
- `_stratified_subsample` fix (BM-04) -- replace `sorted(indices)[:n]` with proportional group downsampling
- Runtime tracking implementation (BM-05) -- add timing around each method in `run_integration_benchmark`
- Parameter provenance storage (BM-06) -- store batch_weight, bio_weight, seed, resolution in DataFrame attrs or as columns
- `_recipe_targeted_panel` scVI fix (BM-07) -- skip normalization when scVI is selected so raw counts are preserved
- Test fixture design (TST-01/02/03) -- synthetic data with controlled properties for unit tests

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope

</user_constraints>

<phase_requirements>

## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| BM-01 | Load pre-computed embeddings from h5ad via h5py (memory <2GB for 2.5M cells) | h5py access pattern for obsm/obs verified; h5ad structure documented |
| BM-02 | Filter NaN rows per-embedding before metric computation | NumPy `np.isnan().any(axis=1)` masking; must align obs arrays |
| BM-03 | Configurable `subsample_n` parameter in `compare_integrations()` | Already exists in `run_full_integration_workflow`; wire to `compare_integrations` |
| BM-04 | Fix `_stratified_subsample` truncation bias | Bug at line 737: `sorted(indices)[:n]` truncates low-index; fix with proportional allocation |
| BM-05 | Add `runtime_s` column to benchmark output | `time.perf_counter()` around method loop in `run_integration_benchmark` |
| BM-06 | Store benchmark parameters alongside results | Use `DataFrame.attrs` (existing pattern, line 424) |
| BM-07 | Fix `_recipe_targeted_panel` -- skip normalization for scVI | Restructure to match `_recipe_visium` pattern (scVI branch before normalize) |
| TST-01 | Unit tests for `compute_integration_metrics` (synthetic, single batch, no celltype) | Existing fixture pattern `_make_batched_adata`; add edge cases |
| TST-02 | Unit tests for `compare_integrations` (NaN embeddings, single method, subsampling) | New test class for NaN/subsample scenarios |
| TST-03 | Unit tests for `_stratified_subsample` (proportionality, n > n_obs, single group) | Extend existing `TestStratifiedSubsample` class |

</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| h5py | 3.15.1 | Direct h5ad file access for obsm/obs | Already installed; avoids full AnnData load |
| numpy | (installed) | NaN masking, array ops | Core dependency |
| pandas | (installed) | DataFrame construction, attrs | Core dependency |
| anndata | (installed) | AnnData objects for test fixtures | Core dependency |
| pytest | >=7.0 | Test framework | Configured in pyproject.toml |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| time (stdlib) | -- | `perf_counter()` for runtime tracking | BM-05 |
| sklearn | (installed) | Metric fallbacks | Already used in integration.py |

No new dependencies required. Everything needed is already installed.

## Architecture Patterns

### Existing Code Structure (DO NOT REORGANIZE)
```
sc_tools/bm/integration.py       # ALL benchmark functions live here
sc_tools/pp/recipes.py            # Preprocessing recipes (_recipe_targeted_panel)
sc_tools/pp/strategy.py           # SmallStrategy (normalization pipeline)
sc_tools/tests/test_integration_benchmark.py  # Existing 32 tests
```

### Pattern 1: h5py Reading from h5ad Files
**What:** Read obsm embeddings and obs categorical columns directly via h5py
**When to use:** BM-01 -- `embedding_files` parameter in `compare_integrations()`
**Example:**
```python
# h5ad stores categoricals as codes + categories groups
import h5py
import numpy as np

def _load_embedding_h5py(path: str, obsm_key: str, obs_keys: list[str]) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Load embedding + obs columns from h5ad without AnnData."""
    with h5py.File(path, "r") as f:
        embedding = f[f"obsm/{obsm_key}"][:]
        obs_data = {}
        for key in obs_keys:
            grp = f[f"obs/{key}"]
            if "categories" in grp:
                # Categorical column
                codes = grp["codes"][:]
                cats = grp["categories"][:]
                if cats.dtype.kind in ("O", "S", "U"):
                    cats = cats.astype(str)
                obs_data[key] = cats[codes]
            else:
                # Plain array
                obs_data[key] = grp[:]
    return embedding, obs_data
```
**Source:** Verified by inspecting actual h5ad file structure (see research step).

**Memory analysis:** A 2.5M-cell dataset with 30-dim embedding = 2.5M x 30 x 4 bytes = 300MB. Plus one obs column ~2.5M x 8 bytes = 20MB. Well under 2GB target.

### Pattern 2: Per-Embedding NaN Masking (BM-02)
**What:** Before computing metrics for each embedding, mask out rows with any NaN
**When to use:** resolVI produces NaN for cells with <5 HVG counts
```python
def _mask_nan_rows(X: np.ndarray, *arrays: np.ndarray) -> tuple[np.ndarray, ...]:
    """Remove rows where embedding has NaN; apply same mask to obs arrays."""
    valid = ~np.isnan(X).any(axis=1)
    n_dropped = (~valid).sum()
    if n_dropped > 0:
        logger.warning("Dropped %d cells with NaN embeddings", n_dropped)
    return (X[valid], *(a[valid] for a in arrays))
```

### Pattern 3: Proportional Stratified Subsampling (BM-04 Fix)
**What:** Current bug: `sorted(indices)[:n]` favors low-index cells. Fix: proportional allocation per group.
**Current buggy code (line 737):**
```python
indices = sorted(indices)[:n]  # BUG: truncates, biases toward low-index groups
```
**Fix pattern:**
```python
# After collecting per-group samples proportionally, the total may slightly
# exceed n due to rounding. Trim randomly instead of by index.
if len(indices) > n:
    indices = rng.choice(indices, size=n, replace=False).tolist()
```

### Pattern 4: DataFrame.attrs for Provenance (BM-06)
**What:** Store benchmark configuration in `DataFrame.attrs` (already used for `scib_fallback`)
**Example:**
```python
df.attrs["benchmark_params"] = {
    "batch_weight": batch_weight,
    "bio_weight": bio_weight,
    "seed": seed,
    "resolution": resolution,
    "scib_backend": use_scib,
    "subsample_n": subsample_n,
}
```
**Note:** `DataFrame.attrs` does not survive all pandas operations (e.g., `pd.concat`). This is acceptable because the benchmark DataFrame is created once and consumed directly.

### Pattern 5: Conditional Normalization in Targeted Panel Recipe (BM-07)
**What:** `_recipe_targeted_panel` currently normalizes before scVI check. Fix: mirror `_recipe_visium` pattern.
**Current buggy flow (recipes.py:244-254):**
```python
normalize_total(adata)  # Always runs -- corrupts raw counts for scVI
log_transform(adata)
strategy.select_features(adata, ...)
if integration == "scvi":
    logger.warning("scVI expects raw counts...")  # Warning but damage is done
```
**Fix pattern (match _recipe_visium):**
```python
if integration == "scvi":
    # scVI needs raw counts -- skip normalize, use seurat_v3 for HVG
    strategy.select_features(adata, n_top_genes=n_top_genes, flavor="seurat_v3", ...)
    ctx = strategy.reduce_and_integrate(adata, integration="scvi", ...)
else:
    normalize_total(adata)
    log_transform(adata)
    strategy.select_features(adata, n_top_genes=n_top_genes, ...)
    ctx = strategy.reduce_and_integrate(adata, integration=integration, ...)
```

### Anti-Patterns to Avoid
- **Loading full AnnData for embeddings:** `ad.read_h5ad()` loads entire X matrix into memory. Use h5py for targeted reads.
- **Silent NaN propagation:** Never pass NaN embeddings to sklearn metrics -- silhouette_score returns NaN, which then propagates through composite scoring.
- **Index-based truncation for subsampling:** `sorted(indices)[:n]` biases toward low-index cells (usually from the first batch/group).

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| h5ad parsing | Custom h5ad reader | h5py with known h5ad structure | h5ad IS HDF5 with a known schema |
| Silhouette score | Custom distance metric | sklearn `silhouette_score` | Already used, well-tested |
| Timing | Custom timer class | `time.perf_counter()` | Standard, nanosecond resolution |
| NaN detection | Manual loop | `np.isnan(X).any(axis=1)` | Vectorized, correct for all dtypes |

## Common Pitfalls

### Pitfall 1: h5ad Categorical Encoding
**What goes wrong:** h5ad stores pandas categoricals as separate `codes` and `categories` datasets. Naively reading `obs/batch` returns a Group, not a Dataset.
**Why it happens:** h5py sees HDF5 groups, not pandas abstractions.
**How to avoid:** Check for `categories` key in group; reconstruct array from `cats[codes]`.
**Warning signs:** TypeError when trying to read obs column as array.

### Pitfall 2: NaN Masking Misalignment
**What goes wrong:** Mask NaN rows from embedding but forget to apply same mask to batch/celltype arrays. Metrics compute on misaligned data.
**Why it happens:** Separate variables, easy to forget one.
**How to avoid:** Single helper function that returns all masked arrays together (Pattern 2 above).
**Warning signs:** ValueError from sklearn about inconsistent sample counts.

### Pitfall 3: Subsample Group Rounding
**What goes wrong:** Proportional allocation `int(len(group) * n / total)` can sum to less than `n` due to floor rounding, or `max(1, ...)` can push total over `n`.
**Why it happens:** Integer division rounding.
**How to avoid:** After proportional allocation, trim or pad to exactly `n` using random selection.
**Warning signs:** Returned subsample size != requested `n`.

### Pitfall 4: Strategy Pattern in _recipe_targeted_panel
**What goes wrong:** `strategy.reduce_and_integrate` with `integration="scvi"` expects raw counts. If `normalize_total` already ran, scVI will train on normalized data and produce poor latent representations.
**Why it happens:** Targeted panel recipe was written before scVI integration was added.
**How to avoid:** Branch on integration method BEFORE normalization (match `_recipe_visium` pattern).
**Warning signs:** scVI training loss unusually high, integration metrics worse than harmony.

### Pitfall 5: DataFrame.attrs Persistence
**What goes wrong:** `DataFrame.attrs` are lost during certain operations (concat, merge, copy in some pandas versions).
**Why it happens:** pandas considers attrs as metadata that may not be preserved.
**How to avoid:** Set attrs after final DataFrame construction. Don't rely on attrs surviving intermediate operations.
**Warning signs:** Missing `benchmark_params` key in downstream consumers.

## Code Examples

### h5ad Structure (Verified)
```
# Actual h5ad internal layout (from research verification):
X                 Dataset  (n_obs, n_vars)  float32
obs/              Group
  _index          Dataset  (n_obs,)         object
  batch/          Group                     # categorical
    categories    Dataset  (n_cats,)        object
    codes         Dataset  (n_obs,)         int8
  celltype/       Group                     # categorical
    categories    Dataset  (n_cats,)        object
    codes         Dataset  (n_obs,)         int8
obsm/             Group
  X_scVI          Dataset  (n_obs, n_latent) float32
```

### Test Fixture Pattern (from existing codebase)
```python
def _make_batched_adata(n_obs=200, n_vars=50, n_batches=3, n_celltypes=4):
    """Existing pattern in test_integration_benchmark.py."""
    rng = np.random.RandomState(42)
    X = rng.randn(n_obs, n_vars).astype(np.float32)
    adata = AnnData(X=X)
    adata.obs["batch"] = [f"batch_{i % n_batches}" for i in range(n_obs)]
    adata.obs["celltype"] = [f"type_{i % n_celltypes}" for i in range(n_obs)]
    adata.obsm["X_good"] = rng.randn(n_obs, 10).astype(np.float32)
    adata.obsm["X_pca"] = rng.randn(n_obs, 10).astype(np.float32)
    return adata
```

### NaN Embedding Fixture (New, for TST-02)
```python
def _make_nan_embedding_adata(n_obs=200, nan_fraction=0.1):
    """AnnData with NaN rows in one embedding (simulates resolVI)."""
    adata = _make_batched_adata(n_obs=n_obs)
    emb = np.random.randn(n_obs, 10).astype(np.float32)
    nan_rows = np.random.choice(n_obs, int(n_obs * nan_fraction), replace=False)
    emb[nan_rows] = np.nan
    adata.obsm["X_resolvi"] = emb
    return adata, nan_rows
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `ad.read_h5ad()` for embeddings | h5py targeted reads | This phase | Memory: full AnnData vs ~300MB for 2.5M cells |
| `celltype_key` only for bio | `bio_key` (any obs column) | This phase | Supports condition, disease_status, etc. |
| `sorted(indices)[:n]` | Proportional group downsampling | This phase | Eliminates truncation bias |

## Open Questions

1. **h5py thread safety for parallel metric computation**
   - What we know: h5py file objects are not thread-safe, but arrays read into memory are fine.
   - What's unclear: Whether parallel metric computation across methods would benefit from threading.
   - Recommendation: Read all embeddings into memory first (they're small), then compute metrics. No threading needed for Phase 1.

2. **bio_key vs celltype_key naming in API**
   - What we know: User decided on `bio_key` parameter that defaults to `celltype_key` value.
   - What's unclear: Whether `bio_key` should replace `celltype_key` in the signature or be an additional parameter.
   - Recommendation: Add `bio_key` as new parameter, deprecate but keep `celltype_key` for backwards compatibility. If both provided, `bio_key` takes precedence.

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | pytest >=7.0 |
| Config file | pyproject.toml `[tool.pytest.ini_options]` |
| Quick run command | `pytest sc_tools/tests/test_integration_benchmark.py -x -v` |
| Full suite command | `pytest sc_tools/tests/ -v --tb=short -q` |

### Phase Requirements -> Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| BM-01 | h5py loading path returns valid embeddings + obs | unit | `pytest sc_tools/tests/test_integration_benchmark.py::TestH5pyLoading -x` | No -- Wave 0 |
| BM-01 | Memory stays <2GB for large embeddings | unit (mock) | `pytest sc_tools/tests/test_integration_benchmark.py::TestH5pyLoading::test_memory_estimate -x` | No -- Wave 0 |
| BM-02 | NaN rows filtered, metrics valid | unit | `pytest sc_tools/tests/test_integration_benchmark.py::TestNaNHandling -x` | No -- Wave 0 |
| BM-03 | subsample_n parameter wired through | unit | `pytest sc_tools/tests/test_integration_benchmark.py::TestCompareIntegrations::test_subsample_n -x` | No -- Wave 0 |
| BM-04 | Proportional group preservation | unit | `pytest sc_tools/tests/test_integration_benchmark.py::TestStratifiedSubsample::test_proportional -x` | No -- Wave 0 |
| BM-05 | runtime_s column in output | unit | `pytest sc_tools/tests/test_integration_benchmark.py::TestRunIntegrationBenchmark::test_runtime_column -x` | No -- Wave 0 |
| BM-06 | benchmark_params in DataFrame.attrs | unit | `pytest sc_tools/tests/test_integration_benchmark.py::TestCompareIntegrations::test_params_in_attrs -x` | No -- Wave 0 |
| BM-07 | scVI path skips normalization | unit | `pytest sc_tools/tests/test_pp.py::TestTargetedPanelScVI -x` | No -- Wave 0 |
| TST-01 | compute_integration_metrics edge cases | unit | `pytest sc_tools/tests/test_integration_benchmark.py::TestComputeMetricsEdgeCases -x` | No -- Wave 0 |
| TST-02 | compare_integrations NaN/single-method | unit | `pytest sc_tools/tests/test_integration_benchmark.py::TestCompareIntegrationsEdgeCases -x` | No -- Wave 0 |
| TST-03 | _stratified_subsample edge cases | unit | `pytest sc_tools/tests/test_integration_benchmark.py::TestStratifiedSubsampleEdgeCases -x` | No -- Wave 0 |

### Sampling Rate
- **Per task commit:** `pytest sc_tools/tests/test_integration_benchmark.py -x -v`
- **Per wave merge:** `pytest sc_tools/tests/ -v --tb=short -q`
- **Phase gate:** Full suite green before `/gsd:verify-work`

### Wave 0 Gaps
- [ ] Test classes for BM-01 through BM-07 and TST-01 through TST-03 (new test cases in existing `test_integration_benchmark.py`)
- [ ] BM-07 tests may go in `test_pp.py` (targeted panel recipe is in `recipes.py`)
- [ ] NaN embedding fixture helper (`_make_nan_embedding_adata`)
- [ ] h5ad fixture helper (write temp h5ad, return path for h5py loading tests)
- No new framework install needed -- pytest already configured

## Sources

### Primary (HIGH confidence)
- Direct code inspection: `sc_tools/bm/integration.py` (all functions, line-by-line)
- Direct code inspection: `sc_tools/pp/recipes.py` (targeted panel recipe)
- Direct code inspection: `sc_tools/pp/strategy.py` (SmallStrategy)
- Direct code inspection: `sc_tools/tests/test_integration_benchmark.py` (32 existing tests)
- h5ad structure verification: created test h5ad and inspected via h5py
- h5py version verified: 3.15.1 installed in current environment

### Secondary (MEDIUM confidence)
- h5ad format specification follows AnnData on-disk format (documented in anndata docs)

### Tertiary (LOW confidence)
- None -- all findings verified against actual code and environment

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- all libraries already installed and verified
- Architecture: HIGH -- patterns derived from direct code inspection
- Pitfalls: HIGH -- bugs identified by reading actual code (line numbers referenced)

**Research date:** 2026-03-20
**Valid until:** 2026-04-20 (stable domain, no external dependencies changing)
