# Phase 8: Multi-Omic Assembly - Research

**Researched:** 2026-03-24
**Domain:** Multi-modal single-cell data integration (MuData, MOFA+, scvi-tools)
**Confidence:** HIGH

## Summary

Phase 8 creates a new `sc_tools/assembly/` module that assembles independently processed per-modality AnnData objects (scRNA, IMC, Visium, Xenium) into a unified MuData object keyed by modality type, with patient-level metadata joins via `subject_id`. The core deliverables are: (1) a `MultiOmicAtlas` class wrapping MuData with hierarchical metadata and cross-modal query methods, (2) pluggable joint embedding using MOFA+ as the recommended default (with MultiVI and TotalVI as optional alternatives), and (3) a new `sct assemble` CLI command group with `build`, `embed`, and `query` subcommands.

The existing codebase provides all necessary integration points: the celltype dispatch pattern (`_base.py` Protocol + registry) serves as the model for joint embedding dispatch, `IOGateway` provides memory-safe loading, `validate_subject_metadata()` provides the join key validation, and the CLI registration pattern from Phase 7's `register_estimate` function shows how to add new command groups.

**Primary recommendation:** Use `mudata>=0.3.1` (scverse standard) as the core container. Implement MOFA+ via `muon.tl.mofa()` as the default joint embedding method. MultiVI (RNA+ATAC only) and TotalVI (RNA+protein only) are niche -- implement them as optional backends but document their modality constraints clearly. MOFA+ is the only method that handles arbitrary modality combinations with missing data, which aligns with the outer join strategy (D-01).

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- **D-01:** Outer join on `subject_id`. Include all patients even with missing modalities -- missing modalities are null in MuData. Filtering to complete cases is easy downstream.
- **D-02:** Hierarchical metadata structure. Patient-level table plus per-modality sample tables. Supports queries at multiple granularities: cell level, sample level (normal vs tumor, per slide), patient level, and potentially organ level.
- **D-03:** N-level abstraction stacking. The stacking is not limited to 3 fixed levels -- supports arbitrary abstraction layers (cell, sample, patient, organ, cohort, etc.) defined by grouping columns in obs metadata.
- **D-04:** Modality-keyed MuData. Keys by modality type: `mdata['rna']`, `mdata['imc']`, `mdata['visium']`, `mdata['xenium']`. Multiple samples per modality are concatenated within each slot. Standard MuData convention.
- **D-05:** Feature spaces kept separate. Each modality keeps its own var (features). Cross-modal analysis happens through shared obs (patient/sample metadata) and joint embeddings in obsm -- not by merging feature spaces.
- **D-06:** Class wrapper -- `MultiOmicAtlas`. Holds MuData + patient/sample metadata tables. Methods like `.patient_view(patient_id)`, `.sample_view(sample_id)`, `.celltype_proportions()`. Encapsulates the multilayer aggregation logic.
- **D-07:** Three methods supported: MOFA+ (factor analysis, handles missing modalities), MultiVI (deep generative, RNA+ATAC), TotalVI (RNA+protein joint model). MOFA+ is the recommended default for multi-platform spatial data.
- **D-08:** Pluggable dispatch pattern. Same architecture as `annotate_celltypes()` -- method registry, runtime selection via `--method` flag. Consistent with existing codebase.
- **D-09:** New `sct assemble` CLI command group with subcommands: `sct assemble build`, `sct assemble embed`, `sct assemble query`.
- **D-10:** V1 query scope: cell type proportions per patient across modalities.

### Claude's Discretion
- How `MultiOmicAtlas` serializes to disk (MuData h5mu format, or separate h5ad + metadata JSON)
- Internal structure of the method registry for joint embedding (enum-based vs class-based dispatch)
- How `sct assemble build` discovers and validates per-modality h5ad files (directory scan vs manifest file vs explicit paths)
- How hierarchical metadata tables are stored in MuData (uns dict, separate DataFrames, or muon's built-in obs mechanism)
- Test fixture design for multi-omic scenarios (synthetic multi-modality data)
- Whether MOFA+, MultiVI, TotalVI wrappers are thin (delegate to library) or thick (handle data prep + postprocessing)

### Deferred Ideas (OUT OF SCOPE)
- Cross-modal marker comparison (e.g., CD20 protein in IMC vs MS4A1 gene in RNA)
- Modality coverage summary table
- Patient-level clinical correlation (survival, treatment response)
- WNN (Weighted Nearest Neighbors) as additional embedding method
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| MOM-01 | Patient/subject metadata join across modalities using validated `subject_id` | MuData outer join via pandas merge on `subject_id` column; reuse `validate_subject_metadata()` for input validation; hierarchical metadata stored in `mdata.uns` as DataFrames |
| MOM-02 | MuData assembly from independently processed per-modality AnnData objects | `mudata.MuData({'rna': adata_rna, 'imc': adata_imc, ...})` constructor; h5mu format for serialization; IOGateway for memory-safe loading of per-modality h5ad files |
| MOM-03 | Joint embedding via multi-modal integration method (MOFA+, MultiVI, TotalVI) | `muon.tl.mofa()` for MOFA+ (default); scvi-tools `MULTIVI`/`TOTALVI` for alternatives; Protocol-based dispatch pattern from celltype module |
| MOM-04 | Cross-modal queries at patient/project level | `MultiOmicAtlas.celltype_proportions()` aggregating across modalities by patient; groupby on configurable abstraction levels |
</phase_requirements>

## Project Constraints (from CLAUDE.md)

- **Orchestrator model:** Default is to dispatch to subagents -- implementation work goes to `pipeline-developer`
- **Output paths:** Active projects live at `~/Documents/projects/active/<project>/` -- never repo root
- **Lazy imports:** Heavy dependencies (scanpy, torch, scvi-tools, muon, mudata) must be loaded at command execution, not at startup. `sct help` must return in <500ms (CLI-08)
- **CLI patterns:** Typer-based, CLIResult envelope, `@cli_handler` decorator, provenance sidecars
- **GPU fallback:** Wrap GPU ops with try/except and fall back to CPU equivalents
- **Linting:** Ruff; never commit failing lint
- **Tests:** Unit + integration; fail-proof with empty/sub/full fixtures

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| mudata | 0.3.3 | MuData container for multi-modal data | scverse standard for multi-modal; h5mu format; AnnData-compatible |
| muon | 0.1.7 | Multi-omics toolkit including MOFA+ wrapper | Official scverse interface for `mofa()`, filtering, processing |
| mofapy2 | 0.7.3 | MOFA+ Python backend | Factor analysis engine; handles missing modalities natively |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| scvi-tools | >=1.0 (already in optional deps) | MultiVI and TotalVI models | When user selects deep generative embedding; already a project dependency |
| h5py | (already installed) | Low-level h5mu/h5ad metadata reading | Pre-load inspection of per-modality files via IOGateway |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| MOFA+ (default) | WNN (Seurat) | WNN is R-only, no Python implementation; deferred to v2 |
| MultiVI | MOFA+ for RNA+ATAC | MOFA+ generalizes; MultiVI is deeper for that specific pair |
| Separate h5ad + JSON | h5mu single file | h5mu is self-contained and standard; recommend h5mu |

**Installation:**
```bash
pip install "mudata>=0.3.1" "muon>=0.1.6" "mofapy2>=0.7.0"
```

**pyproject.toml optional dependency group:**
```toml
# Multi-omic assembly
multiomics = [
    "mudata>=0.3.1",
    "muon>=0.1.6",
    "mofapy2>=0.7.0",
]
```

## Architecture Patterns

### Recommended Module Structure
```
sc_tools/
  assembly/
    __init__.py           # Public API: MultiOmicAtlas, assemble, embed
    _atlas.py             # MultiOmicAtlas class (wraps MuData)
    _build.py             # Assembly logic: load per-modality h5ad, merge obs, create MuData
    _metadata.py          # Hierarchical metadata join (outer join on subject_id, N-level stacking)
    _query.py             # Cross-modal query functions (celltype_proportions, patient_view)
    embed/
      __init__.py         # EmbeddingBackend Protocol + registry
      _base.py            # Protocol definition, register/get functions
      _mofa.py            # MOFA+ backend (muon.tl.mofa wrapper)
      _multivi.py         # MultiVI backend (scvi-tools wrapper)
      _totalvi.py         # TotalVI backend (scvi-tools wrapper)
  cli/
    assemble.py           # sct assemble {build,embed,query} commands
```

### Pattern 1: EmbeddingBackend Protocol (mirrors CelltypeBackend)
**What:** Runtime-checkable Protocol defining the interface all joint embedding backends must satisfy.
**When to use:** Every joint embedding method implements this Protocol.
**Example:**
```python
# Source: existing pattern in sc_tools/tl/celltype/_base.py
from typing import Protocol, runtime_checkable

@runtime_checkable
class EmbeddingBackend(Protocol):
    """Protocol every joint embedding backend must satisfy."""

    @staticmethod
    def run(
        mdata: MuData,
        *,
        n_factors: int = 15,
        **kwargs,
    ) -> tuple[np.ndarray, dict]:
        """
        Run joint embedding.

        Returns
        -------
        embedding
            (n_obs, n_factors) array stored in mdata.obsm.
        metadata
            Dict with method name, parameters, runtime info.
        """
        ...

_BACKENDS: dict[str, type] = {}

def register_embedding_backend(name: str, cls: type) -> None:
    _BACKENDS[name] = cls

def get_embedding_backend(name: str) -> type:
    if name not in _BACKENDS:
        raise ValueError(f"Unknown embedding method '{name}'. Available: {sorted(_BACKENDS)}")
    return _BACKENDS[name]
```

### Pattern 2: MultiOmicAtlas Class
**What:** Wrapper around MuData providing patient-centric access patterns and hierarchical metadata.
**When to use:** The public-facing API for multi-omic analysis.
**Example:**
```python
class MultiOmicAtlas:
    """Multi-omic atlas wrapping MuData with hierarchical metadata."""

    def __init__(self, mdata: MuData, metadata: dict[str, pd.DataFrame] | None = None):
        self.mdata = mdata
        # Hierarchical metadata tables stored in mdata.uns
        self._metadata = metadata or {}

    @classmethod
    def from_modalities(
        cls,
        modalities: dict[str, AnnData],
        subject_key: str = "subject_id",
    ) -> "MultiOmicAtlas":
        """Build atlas from per-modality AnnData objects with outer join on subject_key."""
        mdata = MuData(modalities)
        mdata.update()  # Sync obs across modalities
        # Build hierarchical metadata from obs
        ...
        return cls(mdata, metadata)

    def patient_view(self, patient_id: str) -> MuData:
        """Subset to a single patient across all modalities."""
        mask = self.mdata.obs["subject_id"] == patient_id
        return self.mdata[mask]

    def celltype_proportions(
        self,
        *,
        celltype_key: str = "celltype",
        group_by: str = "subject_id",
    ) -> pd.DataFrame:
        """Cell type proportions per group across all modalities."""
        ...

    def save(self, path: str | Path) -> None:
        """Serialize to h5mu with metadata in uns."""
        self.mdata.write(path)

    @classmethod
    def load(cls, path: str | Path) -> "MultiOmicAtlas":
        """Load from h5mu file."""
        import mudata
        mdata = mudata.read(path)
        return cls(mdata)
```

### Pattern 3: CLI Registration (mirrors register_estimate)
**What:** Register `assemble_app` as a Typer subcommand group.
**When to use:** For the `sct assemble` CLI entry point.
**Example:**
```python
# sc_tools/cli/assemble.py
import typer

assemble_app = typer.Typer(help="Multi-omic assembly commands")

@assemble_app.callback(invoke_without_command=True)
def assemble_callback(ctx: typer.Context) -> None:
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())

@assemble_app.command("build")
@cli_handler
def assemble_build(...):
    """Build MuData from per-modality h5ad files."""
    ...

@assemble_app.command("embed")
@cli_handler
def assemble_embed(...):
    """Run joint embedding on assembled MuData."""
    ...

@assemble_app.command("query")
@cli_handler
def assemble_query(...):
    """Query cross-modal data (e.g., celltype proportions)."""
    ...

# In sc_tools/cli/__init__.py:
# from sc_tools.cli.assemble import assemble_app  # noqa: E402
# app.add_typer(assemble_app, name="assemble")
```

### Pattern 4: MuData Outer Join with Missing Modalities
**What:** How to handle patients with partial modality coverage.
**When to use:** During `sct assemble build`.
**Example:**
```python
# MuData naturally handles different obs sets per modality
# When creating MuData from dict, each modality can have different cells/patients
mdata = MuData({
    "rna": adata_rna,       # Has patients A, B, C
    "imc": adata_imc,       # Has patients A, B only
    "visium": adata_visium,  # Has patients A, C only
})
mdata.update()  # Creates unified obs with NaN for missing modalities

# MOFA+ handles missing modalities natively via use_obs='union'
mu.tl.mofa(mdata, use_obs='union', n_factors=15)
# Result: mdata.obsm['X_mofa'] for all cells, factors learned from available data
```

### Anti-Patterns to Avoid
- **Merging feature spaces:** Never concatenate var across modalities. RNA has genes, IMC has proteins, Visium has spots -- they must remain separate per D-05.
- **Inner join on patients:** D-01 mandates outer join. Never drop patients missing a modality.
- **Loading all modalities at once without memory check:** Use IOGateway for each per-modality h5ad before assembly.
- **Hardcoding abstraction levels:** D-03 requires N-level stacking. Use groupby columns, not fixed hierarchy.
- **Global imports of mudata/muon/mofapy2:** Must be lazy per CLI-08.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Multi-modal data container | Custom dict of AnnData objects | `mudata.MuData` | Handles obs sync, h5mu I/O, modality indexing, scverse ecosystem compatibility |
| MOFA+ factor analysis | Custom matrix factorization | `muon.tl.mofa()` wrapping `mofapy2` | Handles missing data, multi-group, GPU mode, convergence diagnostics |
| MultiVI joint embedding | Custom VAE | `scvi.model.MULTIVI` | Deep generative model with imputation; complex training loop |
| TotalVI joint embedding | Custom CITE-seq VAE | `scvi.model.TOTALVI` | Protein + RNA joint model; library size correction |
| h5mu file format | Custom serialization | `mdata.write()` / `mudata.read()` | HDF5-based, standardized, cross-language support |
| Obs synchronization across modalities | Manual index alignment | `mdata.update()` | Handles union/intersection of observations automatically |

**Key insight:** MuData/muon solve the hard problems (obs alignment, missing modalities, h5mu serialization). The `MultiOmicAtlas` class adds domain-specific convenience methods on top, not reimplementation of container logic.

## Common Pitfalls

### Pitfall 1: MuData.update() Must Be Called After Modifications
**What goes wrong:** Adding modalities or modifying obs without calling `mdata.update()` leaves the global obs/var out of sync with per-modality obs/var.
**Why it happens:** MuData lazily syncs its global annotations with per-modality annotations.
**How to avoid:** Always call `mdata.update()` after adding modalities, subsetting, or modifying obs.
**Warning signs:** `mdata.obs` shape doesn't match sum of per-modality obs shapes.

### Pitfall 2: MOFA+ Requires Normalized Data (Not Raw Counts)
**What goes wrong:** MOFA+ produces poor factors when given raw count matrices.
**Why it happens:** MOFA+ assumes approximately Gaussian data. Raw counts are highly skewed.
**How to avoid:** Each modality should be independently processed (normalized, log-transformed, HVG-selected) before MOFA+ embedding. The phase description says "independently processed" -- enforce this at build time.
**Warning signs:** First factor explains >50% of variance and correlates with library size.

### Pitfall 3: MultiVI/TotalVI Modality Constraints
**What goes wrong:** User tries MultiVI on RNA+IMC or TotalVI on RNA+Visium -- these models are modality-specific.
**Why it happens:** MultiVI only supports RNA+ATAC. TotalVI only supports RNA+protein (CITE-seq style).
**How to avoid:** Validate modality compatibility at dispatch time. MOFA+ is the only method that generalizes across arbitrary modality combinations. Add clear error messages: "MultiVI requires RNA and ATAC modalities. For RNA+spatial, use MOFA+ (default)."
**Warning signs:** Setup errors from scvi-tools when modalities don't match expected types.

### Pitfall 4: Subject ID Inconsistency Across Modalities
**What goes wrong:** Same patient has different `subject_id` values across modalities (e.g., "PAT001" in RNA vs "Patient_001" in IMC).
**Why it happens:** Different labs/platforms use different naming conventions.
**How to avoid:** Validate subject_id consistency at build time using `validate_subject_metadata()`. Require a mapping table if IDs don't match.
**Warning signs:** Outer join produces zero overlap in subject_ids across modalities.

### Pitfall 5: Memory Explosion During Multi-Modality Assembly
**What goes wrong:** Loading 4 large AnnData objects simultaneously exceeds available RAM.
**Why it happens:** Each modality could be hundreds of MB to multiple GB.
**How to avoid:** Use IOGateway to estimate total memory before loading. Load one modality at a time, extract obs metadata first (T1), then load full data (T3) only when needed. Support chunked assembly for very large datasets.
**Warning signs:** System swap usage increases during build.

### Pitfall 6: Celltype Key Inconsistency Across Modalities
**What goes wrong:** Different modalities use different column names for cell types (e.g., "celltype" vs "cell_type" vs "celltype_auto_sctype").
**Why it happens:** Independent processing may use different annotation methods.
**How to avoid:** Require a `celltype_key` parameter in `celltype_proportions()`. Validate that the key exists in each modality's obs before computing.
**Warning signs:** KeyError when computing cross-modal queries.

## Code Examples

### Creating MuData from Per-Modality AnnData
```python
# Source: mudata docs + MuData constructor API
import mudata

# Each adata has subject_id in obs, independently processed
mdata = mudata.MuData({
    "rna": adata_rna,
    "imc": adata_imc,
    "visium": adata_visium,
    "xenium": adata_xenium,
})
mdata.update()  # Sync global obs from per-modality obs

# Access modality
rna = mdata.mod["rna"]  # or mdata["rna"]

# Global obs contains union of all cells
print(mdata.obs.shape)  # (total_cells_across_modalities, ...)

# Save to h5mu
mdata.write("atlas.h5mu")

# Load
mdata = mudata.read("atlas.h5mu")
```

### Running MOFA+ via muon
```python
# Source: muon docs - https://muon.readthedocs.io/en/latest/omics/multi.html
import muon as mu

# Run MOFA+ with union of observations (handles missing modalities)
mu.tl.mofa(
    mdata,
    use_obs="union",       # D-01: outer join, include all cells
    n_factors=15,          # Number of latent factors
    gpu_mode=False,        # Set True if GPU available
)

# Result stored in mdata.obsm["X_mofa"]
embedding = mdata.obsm["X_mofa"]  # (n_total_cells, 15)
```

### N-Level Abstraction Stacking (D-03)
```python
# Hierarchical groupby using arbitrary columns from obs
def aggregate_by_level(
    mdata: MuData,
    *,
    group_cols: list[str],
    value_col: str = "celltype",
) -> pd.DataFrame:
    """Aggregate cell type proportions at any abstraction level."""
    records = []
    for mod_name, mod in mdata.mod.items():
        if value_col not in mod.obs.columns:
            continue
        # Group by requested columns + cell type
        available_cols = [c for c in group_cols if c in mod.obs.columns]
        if not available_cols:
            continue
        counts = mod.obs.groupby(available_cols + [value_col]).size()
        totals = mod.obs.groupby(available_cols).size()
        props = (counts / totals).reset_index(name="proportion")
        props["modality"] = mod_name
        records.append(props)
    return pd.concat(records, ignore_index=True)

# Usage at different levels:
# Cell -> Sample level
aggregate_by_level(mdata, group_cols=["sample_id"], value_col="celltype")

# Cell -> Patient level
aggregate_by_level(mdata, group_cols=["subject_id"], value_col="celltype")

# Cell -> Organ -> Cohort level
aggregate_by_level(mdata, group_cols=["organ", "cohort"], value_col="celltype")
```

### CLI Pattern for sct assemble build
```python
# Source: existing pattern in sc_tools/cli/estimate.py
@assemble_app.command("build")
@cli_handler
def assemble_build(
    inputs: list[str] = typer.Argument(..., help="Per-modality h5ad file paths"),
    modalities: list[str] = typer.Option(None, "--modality", "-m", help="Modality labels (rna, imc, visium, xenium)"),
    output: str = typer.Option("atlas.h5mu", "--output", "-o", help="Output h5mu file path"),
    subject_key: str = typer.Option("subject_id", "--subject-key", help="Column name for subject ID"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate inputs without building"),
    force: bool = typer.Option(False, "--force", help="Bypass memory guard"),
) -> None:
    """Build MuData atlas from per-modality h5ad files."""
    _check_deps(["mudata"])
    from sc_tools.assembly import MultiOmicAtlas
    ...
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Concatenate all modalities into one AnnData | MuData with per-modality slots | mudata 0.1 (2022) | Preserves distinct feature spaces |
| MOFA (R only) | MOFA+ via mofapy2 + muon Python wrapper | 2020-2021 | Native Python support with GPU option |
| Manual obs alignment across modalities | `mdata.update()` auto-sync | mudata 0.2+ | Eliminates manual index juggling |
| scvi-tools setup_anndata for MultiVI | setup_mudata (required) | scvi-tools 1.4+ | MuData is now mandatory input format |
| Custom h5 files for multi-modal | h5mu standard format | mudata 0.1+ | Standardized I/O across scverse ecosystem |

**Deprecated/outdated:**
- `anndata.concat()` for multi-modal: Use MuData instead; concat merges feature spaces which violates D-05
- scvi-tools `setup_anndata` for MultiVI: Replaced by `setup_mudata` in scvi-tools 1.4+
- MOFA2 R package for Python workflows: Use `muon.tl.mofa()` which wraps `mofapy2` natively

## Open Questions

1. **H5mu file size with 4 modalities of large datasets**
   - What we know: h5mu is HDF5-based like h5ad. Each modality stored as a separate group.
   - What's unclear: Total file size when assembling e.g., 500K RNA cells + 100K IMC cells + 200K Visium spots
   - Recommendation: Add a memory/size estimate to `sct assemble build --dry-run` using IOGateway's estimate

2. **MOFA+ convergence with highly disparate modality sizes**
   - What we know: MOFA+ handles missing modalities via `use_obs='union'`. Factor number is configurable.
   - What's unclear: Quality of factors when one modality has 500K cells and another has 10K
   - Recommendation: Log convergence diagnostics (ELBO). Consider subsampling large modalities for factor learning, then projecting.

3. **Hierarchical metadata storage location in MuData**
   - What we know: MuData has `mdata.uns` (dict), `mdata.obs` (global obs), and per-modality `mdata['rna'].obs`
   - What's unclear: Whether patient-level and sample-level DataFrames should go in `uns` or as separate attributes in `MultiOmicAtlas`
   - Recommendation: Store in `mdata.uns['patient_metadata']` and `mdata.uns['sample_metadata']` as DataFrames. This persists through h5mu serialization. The `MultiOmicAtlas` class provides convenience accessors.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| mudata | MOM-02 (MuData container) | No | -- | Must install; no fallback |
| muon | MOM-03 (MOFA+ wrapper) | No | -- | Must install; no fallback |
| mofapy2 | MOM-03 (MOFA+ engine) | No | -- | Must install; no fallback |
| scvi-tools | MOM-03 (MultiVI/TotalVI) | No (locally) | -- | Optional; MOFA+ is default |
| anndata | MOM-02 (per-modality input) | Yes | 0.11.4 | -- |
| scanpy | MOM-01 (preprocessing) | Yes | 1.11.5 | -- |
| pandas | MOM-01/04 (metadata joins, queries) | Yes | (installed) | -- |
| h5py | Pre-load inspection | Yes | (installed) | -- |

**Missing dependencies with no fallback:**
- mudata, muon, mofapy2 -- must be installed as optional dependency group `multiomics`

**Missing dependencies with fallback:**
- scvi-tools (for MultiVI/TotalVI) -- MOFA+ is the default; scvi-tools only needed if user selects `--method multivi` or `--method totalvi`

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | pytest 7.0+ |
| Config file | pyproject.toml `[tool.pytest.ini_options]` |
| Quick run command | `pytest sc_tools/tests/test_assembly.py -x` |
| Full suite command | `pytest sc_tools/tests/ -v --tb=short` |

### Phase Requirements -> Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| MOM-01 | Subject metadata outer join across modalities | unit | `pytest sc_tools/tests/test_assembly.py::test_metadata_outer_join -x` | No -- Wave 0 |
| MOM-01 | Validate subject_id consistency across modalities | unit | `pytest sc_tools/tests/test_assembly.py::test_subject_id_validation -x` | No -- Wave 0 |
| MOM-02 | MuData creation from dict of AnnData objects | unit | `pytest sc_tools/tests/test_assembly.py::test_mudata_build -x` | No -- Wave 0 |
| MOM-02 | MultiOmicAtlas save/load roundtrip (h5mu) | unit | `pytest sc_tools/tests/test_assembly.py::test_atlas_roundtrip -x` | No -- Wave 0 |
| MOM-03 | MOFA+ embedding produces obsm['X_mofa'] | integration | `pytest sc_tools/tests/test_assembly.py::test_mofa_embedding -x` | No -- Wave 0 |
| MOM-03 | Embedding backend dispatch (mofa/multivi/totalvi) | unit | `pytest sc_tools/tests/test_assembly.py::test_embedding_dispatch -x` | No -- Wave 0 |
| MOM-03 | Invalid method name raises ValueError | unit | `pytest sc_tools/tests/test_assembly.py::test_invalid_method -x` | No -- Wave 0 |
| MOM-04 | Cell type proportions per patient across modalities | unit | `pytest sc_tools/tests/test_assembly.py::test_celltype_proportions -x` | No -- Wave 0 |
| MOM-04 | N-level aggregation (sample, patient, organ) | unit | `pytest sc_tools/tests/test_assembly.py::test_nlevel_aggregation -x` | No -- Wave 0 |
| CLI | `sct assemble build` creates h5mu from h5ad inputs | integration | `pytest sc_tools/tests/test_cli_assemble.py::test_build_command -x` | No -- Wave 0 |
| CLI | `sct assemble embed` runs MOFA+ embedding | integration | `pytest sc_tools/tests/test_cli_assemble.py::test_embed_command -x` | No -- Wave 0 |
| CLI | `sct assemble query` returns celltype proportions JSON | integration | `pytest sc_tools/tests/test_cli_assemble.py::test_query_command -x` | No -- Wave 0 |

### Sampling Rate
- **Per task commit:** `pytest sc_tools/tests/test_assembly.py -x`
- **Per wave merge:** `pytest sc_tools/tests/ -v --tb=short`
- **Phase gate:** Full suite green before `/gsd:verify-work`

### Wave 0 Gaps
- [ ] `sc_tools/tests/test_assembly.py` -- covers MOM-01, MOM-02, MOM-03, MOM-04
- [ ] `sc_tools/tests/test_cli_assemble.py` -- covers CLI integration tests for assemble commands
- [ ] Test fixtures: synthetic multi-modality AnnData objects (small: 50 cells/modality, 4 modalities, 3 patients with partial coverage)
- [ ] `mudata` must be added to dev dependencies or test extras for CI

### Test Fixture Design Recommendation
```python
@pytest.fixture
def multi_omic_adatas():
    """4-modality fixture: rna (3 patients), imc (2 patients), visium (2 patients), xenium (1 patient)."""
    rng = np.random.default_rng(42)

    def _make_adata(n_obs, n_vars, patients, modality, var_prefix="GENE"):
        X = rng.negative_binomial(5, 0.3, (n_obs, n_vars)).astype("float32")
        adata = AnnData(X)
        adata.var_names = [f"{var_prefix}{i}" for i in range(n_vars)]
        adata.obs_names = [f"{modality}_cell_{i}" for i in range(n_obs)]
        cells_per_patient = n_obs // len(patients)
        adata.obs["subject_id"] = np.repeat(patients, cells_per_patient)[:n_obs]
        adata.obs["celltype"] = [f"type_{i % 3}" for i in range(n_obs)]
        adata.obs["sample_id"] = [f"{modality}_S{i % 2}" for i in range(n_obs)]
        return adata

    return {
        "rna": _make_adata(60, 200, ["PAT1", "PAT2", "PAT3"], "rna"),
        "imc": _make_adata(40, 40, ["PAT1", "PAT2"], "imc", var_prefix="PROT"),
        "visium": _make_adata(40, 150, ["PAT1", "PAT3"], "visium"),
        "xenium": _make_adata(20, 100, ["PAT1"], "xenium"),
    }
```

## Discretion Recommendations

Based on research, here are recommendations for areas left to Claude's discretion:

### Serialization: Use h5mu
**Recommendation:** Serialize `MultiOmicAtlas` as h5mu via `mdata.write()`. Store hierarchical metadata in `mdata.uns` as DataFrames (pandas DataFrames round-trip through h5mu uns). This is simpler than separate h5ad + JSON and keeps everything in one file.

### Method Registry: Class-Based Dispatch (not enum)
**Recommendation:** Use the same Protocol + dict registry pattern as celltype backends. This is already proven in the codebase, supports third-party extensions via `register_embedding_backend()`, and is more flexible than an enum.

### File Discovery: Explicit Paths (not directory scan)
**Recommendation:** `sct assemble build` takes explicit file paths as positional arguments with `--modality` flags to label them. This is unambiguous, works with any directory layout, and matches the CLI's non-interactive philosophy (CLI-07). Example: `sct assemble build rna.h5ad --modality rna imc.h5ad --modality imc --output atlas.h5mu`.

### Metadata Storage: uns DataFrames
**Recommendation:** Store patient-level and sample-level metadata as `mdata.uns['patient_metadata']` (DataFrame) and `mdata.uns['sample_metadata']` (DataFrame). These survive h5mu round-trips. The `MultiOmicAtlas` class provides `.patient_metadata` and `.sample_metadata` properties that read from uns.

### Wrapper Thickness: Thin Wrappers
**Recommendation:** MOFA+, MultiVI, and TotalVI backends should be thin wrappers: validate inputs, call the library function, extract the embedding, return it in the Protocol's format. Data prep (normalization, HVG selection) should happen before the backend is called, in `_build.py`. This keeps backends swappable and testable.

## Sources

### Primary (HIGH confidence)
- mudata PyPI (0.3.3) - MuData constructor API, h5mu I/O
- [mudata docs](https://mudata.readthedocs.io/en/latest/) - MuData quickstart, API reference
- [muon docs](https://muon.readthedocs.io/en/latest/omics/multi.html) - MOFA+ wrapper (`mu.tl.mofa()`) API and parameters
- [scvi-tools MultiVI docs](https://docs.scvi-tools.org/en/stable/user_guide/models/multivi.html) - MultiVI setup_mudata, training, latent representation
- mofapy2 PyPI (0.7.3) - MOFA+ Python backend version
- Existing codebase: `sc_tools/tl/celltype/_base.py` - Protocol + registry dispatch pattern
- Existing codebase: `sc_tools/cli/__init__.py` - CLI registration pattern, cli_handler decorator
- Existing codebase: `sc_tools/io/gateway.py` - IOGateway tiered loading for memory safety
- Existing codebase: `sc_tools/qc/metadata.py` - validate_subject_metadata() for join key validation

### Secondary (MEDIUM confidence)
- [MOFA+ tutorials](https://biofam.github.io/MOFA2/tutorials.html) - MOFA+ usage patterns, factor interpretation
- [scverse/mudata GitHub](https://github.com/scverse/mudata) - MuData design doc, recent releases
- scvi-tools discourse - setup_mudata requirement for MultiVI in scvi-tools 1.4+

### Tertiary (LOW confidence)
- MOFA+ convergence behavior with highly disparate modality sizes -- based on method description, not empirical testing
- h5mu file size estimates for 4-modality large datasets -- extrapolated from h5ad experience

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - mudata/muon/mofapy2 are scverse standards, versions verified against PyPI
- Architecture: HIGH - mirrors existing proven patterns (celltype dispatch, CLI registration, IOGateway)
- Pitfalls: HIGH - MuData update() requirement, MOFA+ normalization, MultiVI/TotalVI modality constraints are documented in official docs
- MultiVI/TotalVI modality limitations: HIGH - verified from official scvi-tools docs that MultiVI is RNA+ATAC only, TotalVI is RNA+protein only

**Research date:** 2026-03-24
**Valid until:** 2026-04-24 (stable scverse ecosystem, 30-day window)
