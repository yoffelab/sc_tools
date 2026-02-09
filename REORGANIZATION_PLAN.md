# Code Reorganization Plan

## Overview
This document outlines the reorganization strategy to separate generic reusable tools from project-specific analysis code, following scanpy's API pattern.

## Library Name: `sc-tools`

The library follows scanpy's API structure:
- `sc_tools.pl` - Plotting utilities (like `scanpy.pl`)
- `sc_tools.tl` - Analysis tools (like `scanpy.tl`)
- `sc_tools.data` - Data I/O and preprocessing
- `sc_tools.memory` - Memory management and GPU utilities

## New Directory Structure

```
.
├── sc_tools/                        # Generic reusable tools (publication-ready)
│   ├── __init__.py                  # Main package init (imports pl, tl, data, memory)
│   ├── pl/                          # Plotting module (like scanpy.pl)
│   │   ├── __init__.py
│   │   ├── spatial.py               # Spatial plot wrappers
│   │   ├── heatmaps.py              # Heatmap/clustermap utilities
│   │   ├── statistical.py           # Statistical annotations (bars, asterisks)
│   │   └── volcano.py               # Volcano plot utilities
│   ├── tl/                          # Tools module (like scanpy.tl)
│   │   ├── __init__.py
│   │   ├── testing.py               # Mann-Whitney, FDR correction
│   │   ├── colocalization.py        # Correlation, Moran's I, enrichment
│   │   └── deconvolution.py         # Signature gene selection, caching
│   ├── data/                        # Data loading and preprocessing utilities
│   │   ├── __init__.py
│   │   ├── io.py                    # AnnData I/O, caching
│   │   └── preprocessing.py         # QC, normalization, filtering
│   └── memory/                      # Memory management utilities
│       ├── __init__.py
│       ├── profiling.py             # Memory tracking, cleanup
│       └── gpu.py                   # GPU detection and management
│
├── sc_analysis/                     # Project-specific analysis modules
│   ├── __init__.py
│   ├── spatial/                     # Spatial analysis specific to GGO
│   │   ├── __init__.py
│   │   ├── tls.py                   # TLS-specific analysis
│   │   ├── tumor_regions.py         # Tumor region analysis
│   │   └── macrophage.py            # Macrophage localization
│   ├── signatures/                  # Gene signature analysis
│   │   ├── __init__.py
│   │   ├── scoring.py               # Signature scoring
│   │   └── colocalization.py       # Process colocalization
│   └── deconvolution/               # Deconvolution workflows
│       ├── __init__.py
│       └── workflows.py            # Cell2location, Tangram workflows
│
├── scripts/                         # Clean, reproducible analysis scripts
│   ├── phase1_ingestion/           # Phase I: Data ingestion
│   ├── phase2_preprocessing/       # Phase II: QC and integration
│   ├── phase3_deconvolution/        # Phase III: Deconvolution
│   ├── phase4_spatial/              # Phase IV: Spatial analysis
│   └── phase5_visualization/        # Phase V: Visualization
│
├── scripts/dev/                     # Active development/testing scripts
│   ├── test_*.py                   # Testing scripts
│   └── experimental_*.py           # Experimental code
│
└── scripts/archive/                 # Legacy/archived code (read-only)
    ├── old_code/                    # From scripts/old code/
    └── deprecated/                  # Deprecated functions
```

## API Usage Examples

```python
import sc_tools as st

# Plotting (like scanpy.pl)
st.pl.spatial(...)           # Spatial plots
st.pl.heatmaps(...)          # Heatmaps/clustermaps
st.pl.statistical(...)       # Statistical annotations
st.pl.volcano(...)           # Volcano plots

# Tools (like scanpy.tl)
st.tl.mwu(...)               # Mann-Whitney U test
st.tl.colocalization(...)    # Colocalization analysis
st.tl.deconvolution(...)     # Deconvolution utilities

# Data and memory
st.data.io.load_cached_signatures(...)
st.memory.profiling.log_memory(...)
st.memory.gpu.check_gpu_available()
```

## Integration with imc-analysis

The library will integrate functionalities from:
- https://github.com/ElementoLab/imc-analysis
- Relevant IMC analysis tools will be adapted and integrated into `sc_tools`

## Migration Strategy

### Phase 1: Create Structure and Extract Common Code ✅
1. ✅ Create new directory structure
2. ✅ Extract memory profiling utilities → `sc_tools/memory/`
3. ✅ Extract GPU checking → `sc_tools/memory/gpu.py`
4. ✅ Extract caching utilities → `sc_tools/data/io.py`
5. ✅ Extract signature gene selection → `sc_tools/tl/deconvolution.py`

### Phase 2: Extract Plotting Utilities
1. Extract spatial plotting → `sc_tools/pl/spatial.py`
2. Extract heatmap utilities → `sc_tools/pl/heatmaps.py`
3. Extract statistical annotations → `sc_tools/pl/statistical.py`
4. Extract volcano plot code → `sc_tools/pl/volcano.py`

### Phase 3: Extract Statistical Functions
1. Extract Mann-Whitney, FDR → `sc_tools/tl/testing.py`
2. Extract colocalization methods → `sc_tools/tl/colocalization.py`

### Phase 4: Organize Scripts
1. Move test scripts → `scripts/dev/`
2. Move old_code → `scripts/archive/old_code/`
3. Organize production scripts by phase
4. Update imports in all scripts

### Phase 5: Project-Specific Modules
1. Extract TLS analysis → `sc_analysis/spatial/tls.py`
2. Extract macrophage analysis → `sc_analysis/spatial/macrophage.py`
3. Extract signature scoring → `sc_analysis/signatures/scoring.py`

### Phase 6: Integrate imc-analysis Features
1. Review imc-analysis library structure
2. Identify relevant functionalities
3. Adapt and integrate into `sc_tools`

---

## Next immediate steps (plan before executing)

**Goal:** Replicate the same scientific outcomes as before (Phase 3–5 objectives), so the reorganized codebase is **reproducibility-equivalent**. Outputs may differ only in naming (e.g. versioned filenames from `st.pl.save_figure` / `st.data.write_h5ad`).

**1. Ensure Makefile Phase 3–5 works**
- Verify that Makefile targets for Phase III (deconvolution), Phase IV (spatial/niche analysis), and Phase V (visualization) run successfully.
- Confirm that running these phases produces the same results as previously (same figures, same downstream inputs), with versioned filenames where applicable.
- Fix any broken dependencies or paths so that `make phase3`, `make phase4`, `make phase5` (or equivalent) complete without error.

**2. Modular scripts driven by config + sc_tools**
- Refactor analysis scripts so they **do not rewrite shared logic**: use **config files** (e.g. dicts of arguments) and **import functions directly** from `sc_tools` (and later `sc_analysis`).
- Scripts should be thin orchestration: load config → call `st.pl.*` / `st.tl.*` / `st.data.*` with config-driven arguments → write versioned outputs.
- Outcome: easier to navigate, no duplicated implementations; changing behavior happens in one place (sc_tools or config), not inside each script.

Execution order: complete (1) first so the pipeline is reproducible, then proceed to (2) for modularization.

## Benefits

1. **Modularity**: Common code reused across scripts
2. **Maintainability**: Single source of truth for utilities
3. **Testability**: Isolated functions easier to test
4. **Clarity**: Clear separation of generic vs project-specific
5. **Publication-ready**: Clean tools library can be published separately
6. **Development workflow**: Clear dev/archive/production separation
7. **Familiar API**: Follows scanpy's pattern for easy adoption

## Implementation Notes

- All `sc_tools/` code should be well-documented with docstrings
- All `sc_tools/` code should have minimal project-specific dependencies
- Scripts in `scripts/` should import from `sc_tools/` and `sc_analysis/`
- `scripts/dev/` scripts are temporary and should be cleaned up regularly
- `scripts/archive/` is read-only for reference only
- API follows scanpy's pattern: `st.pl.*` for plotting, `st.tl.*` for tools
