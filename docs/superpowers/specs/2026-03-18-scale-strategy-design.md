# Scale Strategy Design Spec v3

**Date:** 2026-03-18
**Status:** Draft v3.1 — post spec-review fixes (correctness bugs, API compatibility, CytoVI)
**Problem:** sc_tools processes all data eagerly in-memory via AnnData. This fails at 7M+ cells (immediate need) and won't scale to atlas-size (10-100M, future).

## Design Summary

Strategy pattern on preprocessing recipes. Two strategies (Small and Large) selected by estimated memory footprint + hardware, overridable in config. LargeStrategy runs cheap ops at full resolution, subsamples for expensive graph + integration operations, projects back via kNN in PCA space. Integration methods are scale-aware within the strategy boundary.

**All platforms in scope.** IMC (dense, ~50 markers, combinatorial panels via MuData) scales to millions/tens of millions of cells and is actually an easier case for LargeStrategy (50 markers × 10M cells = ~2GB dense matrix vs 2000 HVGs × 7M = 56GB).

## 1. Scale Strategy Interface

Phase-based ABC matching the actual recipe structure. Four coordination points where scale decisions matter. State flows between phases via an explicit `SubsampleContext` dataclass, not shared instance attributes.

```python
# sc_tools/pp/strategy.py
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
import numpy as np

@dataclass
class SubsampleContext:
    """Immutable state flowing between strategy phases.
    Created in reduce_and_integrate, consumed in embed_and_cluster."""
    subsample_idx: np.ndarray | None = None
    knn_index: object | None = None  # FAISS index or sklearn BallTree
    use_rep: str = "X_pca"
    subsample_n: int = 1_000_000
    projection_k: int = 30
    random_state: int = 42
    stratify_key: str | None = None  # set from batch_key, not hardcoded


class ScaleStrategy(ABC):
    """Defines how a preprocessing recipe executes based on data scale."""

    name: str

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        pass

    @abstractmethod
    def prepare(self, adata, *, raw_backup: str = "layer",
                filter_patterns: list[str] | None = None) -> None:
        """Backup raw data and filter genes.

        SmallStrategy: adata.raw = adata.copy() (current behavior)
        LargeStrategy: adata.layers["raw_counts"] = adata.X.copy() (avoids 2x memory)
        Both: apply filter_patterns if provided (mito/ribo removal).
        """
        ...

    @abstractmethod
    def select_features(self, adata, *, n_top_genes: int = 2000,
                        batch_key: str | None = None, **kw) -> None:
        """HVG/SVG selection + subsetting. Includes batch_key for multi-sample.
        IMC recipes may skip this (all markers retained)."""
        ...

    @abstractmethod
    def reduce_and_integrate(self, adata, *,
                             integration: str = "harmony",
                             batch_key: str | None = None,
                             n_comps: int = 50, **kw) -> SubsampleContext | None:
        """Normalization + PCA + integration. The scale-critical boundary.

        Returns SubsampleContext if subsampling occurred, None otherwise.
        batch_key is always passed explicitly by the recipe (preprocess()
        defaults to "library_id" and forwards it).

        SmallStrategy: scale -> PCA -> integration on full data
        LargeStrategy:
          - Skip scale() (PCA on log-normalized data avoids densification)
          - Sparse-compatible PCA (arpack/covariance_eigh or randomized SVD)
          - Integration is method-dependent:
            - scVI/CytoVI: full data (mini-batches naturally)
            - Harmony: subsample -> correct -> project corrected PCA back
        """
        ...

    @abstractmethod
    def embed_and_cluster(self, adata, *,
                          ctx: SubsampleContext | None = None,
                          resolution: float = 1.0,
                          n_neighbors: int = 15, **kw) -> None:
        """Neighbors + Leiden + UMAP. Subsample+project in LargeStrategy.

        ctx: SubsampleContext from reduce_and_integrate (reuses kNN index).

        SmallStrategy: neighbors -> leiden -> umap on full data
        LargeStrategy: subsample -> neighbors/leiden/umap -> project back
        """
        ...
```

**Design decisions:**
- **ABC with `@abstractmethod`**, not Protocol. Strategies share real behavior (`__enter__`/`__exit__`), and ABC gives clearer error messages when a method is missing.
- **Phase-based, not function-based.** Four methods map to logical recipe phases. Adding new methods (PHATE, TriMap, scVI embeddings) doesn't change the interface.
- **Integration inside the strategy boundary.** `reduce_and_integrate()` owns both dim reduction and integration because they are entangled at scale (Harmony needs full PCA matrix, scVI mini-batches naturally).
- **`SubsampleContext` dataclass** replaces shared instance state (`self._subsample_idx`, `self._knn_index`). Created by `reduce_and_integrate()`, consumed by `embed_and_cluster()`. Makes data flow explicit, enables testing each phase independently.
- **`batch_key` at strategy level defaults to `None`** — the public `preprocess()` API keeps its existing default of `"library_id"` and forwards it through. No API breakage.
- **Context manager.** `__enter__`/`__exit__` are no-ops for Small/Large but enable DistributedStrategy to manage Dask client lifecycle later — no recipe changes needed.
- **`ResourceSpec` deferred to Phase 2.** Not needed for Phase 1 strategy selection (memory estimate suffices). Avoids speculative complexity.

## 2. Strategy Implementations

### SmallStrategy (default — memory fits comfortably)

Current code path. Delegates to existing `_gpu.py` backend dispatcher (which already auto-detects rapids-singlecell for GPU acceleration). No behavior change from today.

```python
class SmallStrategy(ScaleStrategy):
    name = "small"

    def prepare(self, adata, *, raw_backup="layer", filter_patterns=None):
        adata.raw = adata.copy()  # current behavior, fine at <2M cells
        if filter_patterns:
            filter_genes_by_pattern(adata, patterns=filter_patterns)

    def select_features(self, adata, *, n_top_genes=2000, batch_key=None, **kw):
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes,
                                     batch_key=batch_key, **kw)
        adata._inplace_subset_var(adata.var["highly_variable"])

    def reduce_and_integrate(self, adata, *, integration="harmony",
                              batch_key=None, n_comps=50, **kw):
        scale(adata, max_value=10)
        pca(adata, n_comps=n_comps)       # uses _gpu.py backend
        run_integration(adata, method=integration, batch_key=batch_key, **kw)
        return None  # no subsampling

    def embed_and_cluster(self, adata, *, ctx=None, resolution=1.0,
                          n_neighbors=15, **kw):
        use_rep = _auto_use_rep(adata, use_rep=None)
        neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep)
        leiden(adata, resolution=resolution)
        umap(adata)
```

### LargeStrategy (estimated peak memory exceeds threshold)

Full resolution for cheap ops, subsample for expensive integration + graph ops, project back.

```python
class LargeStrategy(ScaleStrategy):
    name = "large"

    def __init__(self, *, subsample_n=1_000_000, subsample_method="stratified",
                 projection_k=30, random_state=42, has_gpu=False, backed=False):
        self.subsample_n = subsample_n
        self.subsample_method = subsample_method
        self.projection_k = projection_k
        self.random_state = random_state
        self.has_gpu = has_gpu
        self.backed = backed  # True when adata.isbacked — disables .copy() patterns

    def prepare(self, adata, *, raw_backup="layer", filter_patterns=None):
        # Layer-based backup: avoids doubling memory (11GB saved at 7M cells)
        adata.layers["raw_counts"] = adata.X.copy()  # sparse, same footprint
        if filter_patterns:
            filter_genes_by_pattern(adata, patterns=filter_patterns)

    def select_features(self, adata, *, n_top_genes=2000, batch_key=None, **kw):
        # HVG on full sparse matrix — expensive but feasible (sparse stays sparse)
        # IMC recipes skip this entirely (all markers retained)
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes,
                                     batch_key=batch_key, **kw)
        adata._inplace_subset_var(adata.var["highly_variable"])

    def reduce_and_integrate(self, adata, *, integration="harmony",
                              batch_key=None, n_comps=50, **kw) -> SubsampleContext | None:
        # SKIP scale() — avoids densification (7M x 2000 x 4B = 56GB)
        # Sparse-compatible PCA
        _sparse_pca(adata, n_comps=n_comps, gpu=self.has_gpu)

        # Integration is method-dependent at scale:
        if integration in ("scvi", "cytovi"):
            # scVI/CytoVI mini-batch naturally — run on full data
            if integration == "cytovi":
                run_cytovi(adata, batch_key=batch_key, **kw)
            else:
                run_scvi(adata, batch_key=batch_key, **kw)
            # Return ctx with use_rep so embed_and_cluster uses latent space
            use_rep = _auto_use_rep(adata, use_rep=None)  # detects X_scVI/X_cytovi
            return SubsampleContext(
                subsample_n=self.subsample_n,
                projection_k=self.projection_k,
                random_state=self.random_state,
                stratify_key=batch_key,
                use_rep=use_rep,
            )  # subsample_idx=None signals "no subsampling done yet"

        # Harmony needs full matrix — subsample, correct, project back
        ctx = SubsampleContext(
            subsample_n=self.subsample_n,
            projection_k=self.projection_k,
            random_state=self.random_state,
            stratify_key=batch_key,
        )
        ctx.subsample_idx = subsample_stratified(
            adata, n=ctx.subsample_n, stratify_key=batch_key,
            random_state=ctx.random_state
        )
        adata_sub = adata[ctx.subsample_idx].copy()
        run_integration(adata_sub, method=integration,
                        batch_key=batch_key, **kw)
        # Project corrected representation back to full data
        ctx.use_rep = _auto_use_rep(adata_sub)
        ctx.knn_index = build_knn_index(adata_sub, use_rep=ctx.use_rep)
        project_representation(adata, ctx.knn_index, adata_sub,
                               rep_key=ctx.use_rep, k=ctx.projection_k)
        return ctx

    def embed_and_cluster(self, adata, *, ctx=None, resolution=1.0,
                          n_neighbors=15, **kw):
        # Build SubsampleContext if not provided (standalone call)
        if ctx is None:
            batch_key = kw.get("batch_key")
            ctx = SubsampleContext(
                subsample_n=self.subsample_n,
                projection_k=self.projection_k,
                random_state=self.random_state,
                stratify_key=batch_key,
                use_rep=_auto_use_rep(adata, use_rep=None),
            )

        # Subsample if not already done (scVI path has ctx but no subsample_idx)
        if ctx.subsample_idx is None:
            ctx.subsample_idx = subsample_stratified(
                adata, n=ctx.subsample_n, stratify_key=ctx.stratify_key,
                random_state=ctx.random_state
            )

        adata_sub = adata[ctx.subsample_idx].copy()
        # Use the correct representation (X_pca, X_scVI, X_cytovi, etc.)
        neighbors(adata_sub, n_neighbors=n_neighbors, use_rep=ctx.use_rep)
        leiden(adata_sub, resolution=resolution)
        umap(adata_sub)

        # Build kNN index if not reused from reduce_and_integrate
        if ctx.knn_index is None:
            ctx.knn_index = build_knn_index(adata_sub, use_rep=ctx.use_rep)

        project_labels(adata, ctx.knn_index, adata_sub,
                       label_key="leiden", k=ctx.projection_k)
        project_umap(adata, ctx.knn_index, adata_sub,
                     k=ctx.projection_k)

        # Mark representative cells and store metadata
        adata.obs["_is_representative"] = False
        adata.obs.loc[adata.obs_names[ctx.subsample_idx],
                      "_is_representative"] = True
        _store_projection_info(adata, ctx)
```

**Why no MediumStrategy:** The existing `_gpu.py` backend dispatcher already auto-detects rapids-singlecell. SmallStrategy with GPU is functionally identical to what MediumStrategy would do. Eliminating it reduces complexity with no capability loss.

**Phase 1 integration scope:** Harmony, scVI, and CytoVI. CytoVI (existing IMC integration) mini-batches like scVI — runs on full data. BBKNN and ComBat are excluded — BBKNN's batch-balanced graph construction doesn't have a clear subsample+project path, and ComBat's linear model assumptions break with subsampling. These can be added in Phase 2 after empirical validation.

## 3. Strategy Selection & Configuration

### Memory-based auto-selection

Thresholds are based on estimated peak memory, not cell count alone. 7M × 300 genes (Xenium) is very different from 7M × 2000 HVGs (Visium HD). Platform-specific density tables inform the estimate.

```python
from sc_tools.memory.profiling import estimate_adata_memory

# Platform-specific density estimates for memory calculation
PLATFORM_DENSITY = {
    # platform: (typical n_vars after HVG, sparsity, dtype_bytes)
    "visium_hd": (2000, 0.85, 4),   # RNA, sparse, float32
    "xenium": (300, 0.70, 4),       # targeted panel, moderately sparse
    "cosmx": (1000, 0.80, 4),       # targeted panel, sparse
    "imc": (50, 0.0, 4),            # dense protein, all markers retained
    "merfish": (500, 0.75, 4),      # targeted panel
}

def select_strategy(adata, config=None, platform=None) -> ScaleStrategy:
    """Pick strategy based on estimated memory footprint + hardware + config."""

    # Config override takes priority
    if config and config.get("backend") not in (None, "auto"):
        return _build_strategy(config["backend"], config)

    # Use platform density table if available, else estimate from data
    if platform and platform in PLATFORM_DENSITY:
        n_vars_est, sparsity, dtype_bytes = PLATFORM_DENSITY[platform]
    else:
        n_vars_est = config.get("n_top_genes", 2000) if config else 2000
        sparsity = _estimate_sparsity(adata.X)
        dtype_bytes = 4

    # Peak memory: dense worst case for scale() in SmallStrategy
    dense_peak_gb = (adata.n_obs * n_vars_est * dtype_bytes) / 1e9
    available_gb = _get_available_memory_gb()  # system or GPU
    has_gpu = check_gpu_available() and _has_rapids()

    # If dense matrix would exceed 50% of available memory → LargeStrategy
    if dense_peak_gb > available_gb * 0.5:
        log.info(f"Estimated peak {dense_peak_gb:.1f}GB exceeds 50% of "
                 f"{available_gb:.0f}GB available → LargeStrategy")
        return LargeStrategy(has_gpu=has_gpu, **_get_large_config(config))

    # Backed AnnData always routes to Large (cannot .copy() freely)
    if adata.isbacked:
        return LargeStrategy(has_gpu=has_gpu, backed=True,
                             **_get_large_config(config))

    return SmallStrategy()
```

### Config integration

```yaml
# project config.yaml or recipe kwargs
preprocessing:
  backend: auto          # auto | small | large
  platform: xenium       # optional — inferred from modality if not set
  large_strategy:
    subsample_n: 1000000
    subsample_method: stratified  # stratified | random | leverage
    projection_k: 30
    random_state: 42
```

- `auto` is default — inspects memory footprint + hardware
- `backend: large` forces subsample+project (useful for testing the pipeline before running on 7M)
- `backend: small` forces current behavior (escape hatch)
- `subsample_method: leverage` — leverage-score sampling (Seurat v5 approach) as future option; `stratified` is default for Phase 1

### Public API — `preprocess()` stays stable

The public `preprocess()` signature is **unchanged**. Strategy selection is internal. The only new parameter is an optional `config` dict forwarded from project config.

```python
# Existing signature — preserved exactly
def preprocess(
    adata: AnnData,
    modality: str = "visium",
    batch_key: str = "library_id",  # default preserved, no regression
    integration: str = "scvi",
    ...
    **kwargs,
) -> AnnData:
    # NEW: select strategy based on data size + config
    strategy = select_strategy(adata, config=kwargs.pop("config", None),
                               platform=modality)
    # Dispatch to recipe with strategy
    recipe_fn = RECIPE_DISPATCH[modality]
    return recipe_fn(adata, batch_key=batch_key, integration=integration,
                     strategy=strategy, ...)
```

### Recipe integration

Recipe signatures gain a `strategy` parameter. All other parameters stay the same — no breaking changes.

```python
def _recipe_visium(adata, batch_key, integration, n_top_genes, resolution,
                   filter_patterns, use_gpu, strategy=None, **kwargs):
    """Existing signature + strategy parameter."""
    if strategy is None:
        strategy = SmallStrategy()  # backward compat for direct calls
    with strategy:
        strategy.prepare(adata, filter_patterns=filter_patterns)
        strategy.select_features(adata, n_top_genes=n_top_genes,
                                  flavor="seurat_v3", batch_key=batch_key)
        ctx = strategy.reduce_and_integrate(
            adata, integration=integration, batch_key=batch_key,
            n_comps=kwargs.get("n_comps", 50))
        strategy.embed_and_cluster(adata, ctx=ctx,
                                    resolution=resolution)

def _recipe_imc(adata, batch_key, integration, n_top_genes, resolution,
                filter_patterns, use_gpu, strategy=None, **kwargs):
    """IMC recipe — dense protein data, pluggable normalization."""
    if strategy is None:
        strategy = SmallStrategy()
    with strategy:
        strategy.prepare(adata, filter_patterns=filter_patterns)
        # IMC: no HVG selection — all markers retained
        # Normalization is pluggable (arcsinh, CLR, log1p, etc.)
        normalize_imc(adata, method=kwargs.get("normalization", "arcsinh"))
        ctx = strategy.reduce_and_integrate(
            adata, integration=integration, batch_key=batch_key,
            n_comps=min(kwargs.get("n_comps", 50), adata.n_vars - 1))
        strategy.embed_and_cluster(adata, ctx=ctx,
                                    resolution=resolution)
```

**IMC note:** With ~50 markers and dense matrix, IMC is actually the easiest case for LargeStrategy. 50 markers × 10M cells = ~2GB dense — well within GPU memory. The memory threshold handles this automatically: IMC only routes to LargeStrategy when cell count is genuinely large.

**Pluggable normalization:** IMC normalization method is a recipe parameter (`normalization="arcsinh"` default), not hardcoded. The recipe owns the normalization call, not the strategy. This allows experimenting with CLR, log1p, or other methods without touching strategy code.

Recipe doesn't know whether it's processing 50K or 7M cells. Strategy handles memory management, subsampling, projection, and integration scaling.

## 4. Projection Module

```python
# sc_tools/pp/projection.py

def subsample_stratified(adata, n=1_000_000, stratify_key=None,
                         random_state=42) -> np.ndarray:
    """Proportional stratified subsample with minimum group representation.
    Falls back to random if stratify_key is None or missing from obs.
    Handles zero-count categorical groups gracefully (skipped, not errored).
    If n >= adata.n_obs, returns all indices (no subsampling)."""

def build_knn_index(adata_sub, use_rep="X_pca"):
    """Build FAISS index on subsample representation.
    GPU if available, CPU fallback. Returns reusable index."""

def project_labels(adata, knn_index, adata_sub,
                   label_key="leiden", k=30):
    """kNN majority vote label transfer.
    Stores confidence in adata.obs[f"{label_key}_confidence"]."""

def project_umap(adata, knn_index, adata_sub, k=30):
    """Weighted average UMAP coordinate transfer.
    Reuses kNN index from project_labels."""

def project_representation(adata, knn_index, adata_sub,
                           rep_key="X_pca_harmony", k=30):
    """Project corrected embedding back to full data.
    Used by LargeStrategy.reduce_and_integrate for Harmony."""
```

**Implementation details:**
- kNN: FAISS-GPU if available → FAISS-CPU → sklearn BallTree (graceful fallback chain)
- FAISS is a **soft dependency** — imported at call time with clear error message
- Index built once, reused for all projections (labels, UMAP, corrected embeddings)
- FAISS-GPU: 7M queries against 1M index ≈ 10-30s. CPU fallback: 5-15min with batching.

**Reproducibility:**
```python
# Subsample membership stored in obs["_is_representative"] (boolean mask).
# uns stores parameters only — no 1M-element lists.
adata.uns["_projection_info"] = {
    "strategy": "large",
    "subsample_n": 1_000_000,
    "subsample_method": "stratified",
    "stratify_key": "library_id",  # actual key used, not hardcoded
    "projection_k": 30,
    "projection_method": "knn_faiss",
    "use_rep": "X_pca_harmony",    # representation used for kNN
    "random_state": 42,
    "median_confidence": 0.87,
    "low_confidence_pct": 0.03,
    "sc_tools_version": sc_tools.__version__,
    "timestamp": datetime.now(timezone.utc).isoformat(),
}
```

**QC guardrails:**
- If `low_confidence_pct > 0.10`, warn suggesting larger subsample or checking for underrepresented populations
- Confidence score per cell enables downstream filtering

## 5. Sparse-Compatible PCA

LargeStrategy must avoid densifying the expression matrix. Two paths:

```python
def _sparse_pca(adata, n_comps=50, gpu=False):
    """PCA that operates on sparse matrices without densification.

    Priority:
    1. GPU + rapids: cuml.decomposition.TruncatedSVD on CuPy sparse
       (VERIFY: cuML TruncatedSVD sparse support — needs empirical test)
    2. scanpy with arpack/covariance_eigh solver (handles zero-centering on sparse)
    3. sklearn.decomposition.TruncatedSVD(algorithm="randomized") on scipy sparse
    """
    if gpu and _has_rapids():
        # cuML path — needs verification that sparse input is supported
        try:
            _cuml_truncated_svd(adata, n_comps=n_comps)
            return
        except TypeError:
            log.warning("cuML TruncatedSVD doesn't support sparse input; "
                        "falling back to scanpy arpack solver")

    # scanpy's PCA with sparse-compatible solver
    # arpack and covariance_eigh handle zero-centering without densification
    sc.pp.pca(adata, n_comps=n_comps, svd_solver="arpack", zero_center=True)
```

**Key insight:** scanpy's `arpack` and `covariance_eigh` solvers already handle zero-centering on sparse matrices without densification. This is simpler and better tested than sklearn's `TruncatedSVD` (which skips centering). The GPU path via cuML `TruncatedSVD` needs empirical verification for sparse input support before relying on it.

## 6. Memory Management at 7M Cells

Critical memory analysis that drives LargeStrategy design:

| Operation | SmallStrategy | LargeStrategy | Savings |
|-----------|--------------|---------------|---------|
| Raw backup | `adata.raw = adata.copy()` (+11GB) | `adata.layers["raw_counts"]` (sparse ref) | ~11GB |
| `scale()` | Densifies: 7M × 2000 × 4B = **56GB** | **Skipped** — PCA on log-normalized sparse | ~56GB |
| PCA | Full SVD on dense matrix | arpack on sparse (or GPU TruncatedSVD) | ~40GB working mem |
| Harmony | Full 7M × 50 PCA matrix | Subsample 1M → correct → project back | ~1GB vs ~2.6GB |
| Neighbors | Full 7M-cell graph | 1M-cell graph + kNN projection | Orders of magnitude |

**IMC at scale (10M cells × 50 markers):**

| Operation | SmallStrategy | LargeStrategy | Notes |
|-----------|--------------|---------------|-------|
| Raw backup | +2GB | layers: same footprint | Dense, but small per-cell |
| Normalize | arcsinh/CLR on 2GB dense | Same (already dense, fast) | No densification penalty |
| PCA | Full SVD on 2GB | arpack on 2GB | Feasible either way |
| Neighbors | 10M-cell graph: **expensive** | 1M-cell graph + project | Main benefit for IMC |

IMC's bottleneck at scale is graph construction and clustering, not memory. LargeStrategy helps via subsampling the graph operations.

## 7. File Layout

### New files
| File | Purpose | ~Lines |
|------|---------|--------|
| `sc_tools/pp/strategy.py` | ABC + SubsampleContext + 2 implementations + selector | ~400 |
| `sc_tools/pp/projection.py` | Subsample + kNN index + project functions | ~250 |

### Modified files
| File | Change | Scope |
|------|--------|-------|
| `sc_tools/pp/recipes.py` | Wrap each recipe in `with strategy:`, delegate to 4 phase methods, pass `ctx` between phases. 4 recipe functions (visium, xenium, cosmx, imc); `visium_hd_cell` dispatches to xenium recipe. | ~50-80 lines per recipe |
| `sc_tools/pp/normalize.py` | Add `normalize_imc(method=...)` dispatcher for pluggable IMC normalization | ~30 lines |

### Untouched files
`_gpu.py`, `reduce.py`, `integrate.py`, `pipeline.py`, `memory/gpu.py` — strategies call these internally.

### Dependencies
```toml
# pyproject.toml
[project.optional-dependencies]
scale = ["faiss-cpu>=1.7"]       # kNN projection
scale-gpu = ["faiss-gpu>=1.7"]   # GPU-accelerated kNN
# rapids-singlecell already in [gpu] extra
```

All new deps are soft — imported at call time with clear error messages.

## 8. Helper Functions (to implement)

Functions referenced in strategy/projection code that do not yet exist:

| Function | Location | Purpose |
|----------|----------|---------|
| `_has_rapids()` | `strategy.py` | Returns `True` if `rapids_singlecell` is importable. Wraps existing try/except pattern from `_gpu.py`. |
| `_get_available_memory_gb()` | `strategy.py` | Returns available system RAM (or GPU VRAM if GPU path). Uses `psutil.virtual_memory()` / `torch.cuda.mem_get_info()`. |
| `_estimate_sparsity(X)` | `strategy.py` | Returns fraction of zeros in sparse/dense matrix. Fallback when platform not in density table. |
| `_build_strategy(name, config)` | `strategy.py` | Factory: maps `"small"` → `SmallStrategy()`, `"large"` → `LargeStrategy(**config)`. |
| `_get_large_config(config)` | `strategy.py` | Extracts `subsample_n`, `projection_k`, etc. from config dict with defaults. |
| `_store_projection_info(adata, ctx)` | `strategy.py` | Writes `adata.uns["_projection_info"]` and `adata.obs["_is_representative"]` from `SubsampleContext`. |
| `normalize_imc(adata, method)` | `normalize.py` | Dispatcher: `"arcsinh"` → `arcsinh_transform()`, `"clr"` → CLR, etc. Pluggable IMC normalization. |

Existing functions used (no changes needed):
- `_auto_use_rep(adata, use_rep)` — `reduce.py:33`, auto-detects `X_scVI`/`X_pca_harmony`/`X_pca`
- `check_gpu_available()` — `memory/gpu.py`
- `estimate_adata_memory()` — `memory/profiling.py:149`

## 9. Testing Strategy

### Unit tests (`test_strategy.py`)
- Strategy selection by memory estimate + hardware + config override
- `select_strategy()` routes backed AnnData to LargeStrategy
- Platform density table influences threshold correctly (IMC 50 vars vs Visium 2000 vars)
- SmallStrategy `prepare()` uses `.copy()`, LargeStrategy uses layers
- LargeStrategy `reduce_and_integrate()` subsamples for Harmony but not scVI
- `reduce_and_integrate()` returns `SubsampleContext` for Harmony, `None` for scVI
- `SubsampleContext` flows correctly from `reduce_and_integrate` to `embed_and_cluster`
- Config override forces strategy regardless of data size
- Boundary cases: just-below-threshold, no GPU available

### Unit tests (`test_projection.py`)
- Stratified subsample proportions and minimum representation per group
- Zero-count categorical groups in stratify_key are skipped gracefully
- `random_state` produces deterministic subsamples
- Label projection on well-separated clusters (100% confidence)
- Label projection on overlapping clusters (lower confidence, verify range [0,1])
- UMAP projection continuity (within convex hull of neighbors)
- kNN index built once, reused across project_labels and project_umap
- `_projection_info` contains all required keys; `_is_representative` mask in obs matches subsample
- Soft dependency error message when FAISS not installed

### Synthetic edge cases
- Zero-count categorical library_id (exists in obs categories but 0 cells)
- Single-sample dataset (no batch dimension for stratification)
- Dataset just above/below memory threshold
- Empty subsample group after filtering
- IMC-shaped data (50 vars, dense, 100K cells) routes correctly
- `subsample_n > adata.n_obs` → returns all indices (no subsampling)
- scVI path: `embed_and_cluster` uses `X_scVI` for neighbors, not `X_pca`

### Integration test (real data)
- Run `_recipe_visium` with `backend: large` on subsampled real checkpoint (~50K cells, stratify to ~10K)
- Verify all outputs: `X_pca`, `X_umap`, `leiden`, `leiden_confidence`, `_is_representative`, `_projection_info`, `layers["raw_counts"]`
- Verify `scale()` was NOT called (no dense matrix in memory profile)
- Verify scVI integration runs on full data, Harmony on subsample

### Validation experiment (required before production use)
- Take a ~500K cell subset of the 7M spatial dataset
- Run both SmallStrategy (ground truth) and LargeStrategy on identical data
- Compare: ARI between leiden clusterings, Procrustes distance on UMAP, per-cell label concordance
- Target: >95% label concordance (Seurat v5 / Symphony report >97%)
- Document results in project Journal.md

### Harmony quality validation plan
- Harmony subsample+project is less validated than scVI at scale
- Run Harmony on 1M subsample, project to full 7M, compare iLISI/cLISI metrics
- If batch correction quality degrades significantly, recommend scVI for this dataset
- Document threshold: if iLISI drops >20% vs full-data Harmony on 500K subset, flag

## 10. HPC Deployment (cayuga)

SLURM job sizing recommendations for the 7M spatial dataset:

```bash
# Phase 1: LargeStrategy with Harmony
#SBATCH --partition=scu-gpu
#SBATCH --gres=gpu:a100:1        # 80GB VRAM for rapids PCA + FAISS
#SBATCH --mem=128G               # 7M cells × 300 genes sparse + overhead
#SBATCH --cpus-per-task=16
#SBATCH --time=4:00:00

# Phase 1: LargeStrategy with scVI (full data)
#SBATCH --partition=scu-gpu
#SBATCH --gres=gpu:a100:1        # scVI needs GPU memory for mini-batches
#SBATCH --mem=192G               # scVI + full data in memory
#SBATCH --cpus-per-task=16
#SBATCH --time=8:00:00
```

## 11. Migration Path

### Phase 1 (immediate — 7M spatial dataset)
- Implement `strategy.py` (ABC + SubsampleContext + SmallStrategy + LargeStrategy)
- Implement `projection.py` (subsample, kNN, project functions)
- Wire into recipes (4 recipes × ~60 lines each)
- Add pluggable IMC normalization to `normalize.py`
- Integration scope: Harmony, scVI, CytoVI (BBKNN/ComBat deferred)
- Test on 7M spatial dataset on cayuga (A100)
- Run validation experiment (SmallStrategy vs LargeStrategy on 500K subset)

### Phase 2 (near-term)
- `ResourceSpec` for HPC job auto-sizing
- BBKNN + ComBat LargeStrategy support (after empirical validation)
- Leverage-score sampling option (Seurat v5 approach)
- Geometric subsampling (spatially-aware — preserves spatial density)
- Parametric UMAP option (train on subsample, transform full dataset)
- Projection QC metrics in HTML reports
- RMM pool configuration in `memory/gpu.py` for rapids stability

### Phase 3 (future — atlas scale)
- `DistributedStrategy` composing with LargeStrategy:
  - Zarr v3 checkpoints → `read_lazy()` → Dask-backed AnnData
  - Cheap ops: scanpy on Dask arrays or rapids on Dask-CuPy arrays
  - Expensive ops: same subsample+project (subsample size is fixed regardless of total)
  - Dask client managed via `__enter__`/`__exit__` context manager
- DuckDB on Parquet obs for cross-project meta-analysis
- Incremental algorithms (mean, variance, HVG) à la CZI Census for datasets exceeding single-node memory
- Zarr v3 + Parquet checkpoint format (separate spec)

**Note:** Zarr v3 + Parquet obs export scoped to a separate spec. Keeping this spec focused on the strategy pattern and projection module.

## 12. References

Architecture informed by:
- NVIDIA rapids-singlecell: GPU-accelerated single-cell at scale (arXiv 2603.02402)
- Tahoe-100M: 100M cell foundation model atlas (bioRxiv 2025.02.20.639398)
- CZI CELLxGENE Census: TileDB-SOMA for 65M+ cells (NAR 2025)
- LatchBio: anntiles + PMTiles + DuckDB for spatial atlas visualization
- AnnData 0.12+: read_lazy(), Zarr v3, Dask-backed arrays
- scanpy 1.11+: experimental Dask support for >100M cells
- scprocess: atlas-scale Snakemake pipeline (bioRxiv 2026.03.09)
- Seurat v5 / BPCells: sketch-based analysis with leverage-score sampling
- Symphony: reference mapping with >97% label transfer accuracy
