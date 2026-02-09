# Reorganization Status

## ✅ Completed

### 1. Directory Structure Created
- ✅ `sc_tools/` - Generic reusable tools (renamed from ggo_tools)
  - ✅ `pl/` - Plotting module (like scanpy.pl)
    - ✅ `spatial.py` - plot_spatial_plain_he, plot_spatial_categorical, plot_spatial_continuous
    - ✅ `heatmaps.py` - hex_to_rgb, cluster_within_groups, annotation_colors_from_categories, signature_score_heatmap
    - ✅ `statistical.py` - get_asterisk, add_significance_bars, plot_boxplot_with_stats
    - ✅ `volcano.py` - volcano_plot, volcano_plot_faceted
    - ✅ `save.py` - save_figure (PDF+PNG, versioned, dpi=300)
  - ✅ `tl/` - Tools module (like scanpy.tl)
    - ✅ `testing.py` - Mann-Whitney U, FDR correction
    - ✅ `colocalization.py` - pearson_correlation, morans_i, morans_i_batch, neighborhood_enrichment (squidpy)
    - ✅ `deconvolution.py` - Signature gene selection with caching
  - ✅ `data/` - Data I/O utilities
    - ✅ `io.py` - Caching utilities
  - ✅ `memory/` - Memory management utilities
    - ✅ `profiling.py` - Memory tracking, cleanup
    - ✅ `gpu.py` - GPU detection and management
- ✅ `sc_analysis/` - Project-specific analysis (renamed from ggo_analysis)
- ✅ `scripts/dev/` - Development/testing scripts
- ✅ `scripts/archive/` - Legacy code
- ✅ `scripts/phase*/` - Organized by analysis phase

### 2. API Structure
- ✅ Main package `__init__.py` imports `pl`, `tl`, `data`, `memory`
- ✅ `sc_tools.pl` module structure (like scanpy.pl)
- ✅ `sc_tools.tl` module structure (like scanpy.tl)
- ✅ Basic functions in `tl.testing`: `mwu()`, `fdr_correction()`

### 3. Plotting Utilities (`sc_tools/pl/`)
- ✅ `spatial.py` - plot_spatial_plain_he, plot_spatial_categorical, plot_spatial_continuous
- ✅ `heatmaps.py` - hex_to_rgb, cluster_within_groups, annotation_colors_from_categories, signature_score_heatmap
- ✅ `statistical.py` - get_asterisk, add_significance_bars, plot_boxplot_with_stats
- ✅ `volcano.py` - volcano_plot, volcano_plot_faceted
- ✅ `save.py` - save_figure (PDF+PNG, versioned, dpi=300)

### 4. TL Functions (`sc_tools/tl/`)
- ✅ `testing.py` - Mann-Whitney U, FDR correction
- ✅ `colocalization.py` - pearson_correlation, morans_i, morans_i_batch, neighborhood_enrichment (squidpy)
- ✅ `io` - write_h5ad (versioned by default)

## 🚧 In Progress

- (Optional) Refactor existing scripts to use `st.pl.spatial`, `st.pl.heatmaps`, `st.pl.save_figure`

## 📋 To Do

### 5. Project-Specific Modules (`sc_analysis/`)
- [ ] `spatial/tls.py` - TLS-specific analysis
- [ ] `spatial/tumor_regions.py` - Tumor region analysis
- [ ] `spatial/macrophage.py` - Macrophage localization
- [ ] `signatures/scoring.py` - Signature scoring
- [ ] `signatures/colocalization.py` - Process colocalization
- [ ] `deconvolution/workflows.py` - Cell2location, Tangram workflows

### 6. Script Organization
- [ ] Move test scripts → `scripts/dev/`
- [ ] Move `scripts/old_code/` → `scripts/archive/old_code/`
- [ ] Organize production scripts by phase
- [ ] Update imports in all scripts to use `sc_tools.pl.*`, `sc_tools.tl.*`

### 7. Integration with imc-analysis
- [ ] Review imc-analysis library structure
- [ ] Identify relevant functionalities
- [ ] Adapt and integrate into `sc_tools`

### 8. Documentation
- [ ] Update `Architecture.md` with new structure
- [ ] Create migration guide for updating existing scripts
- [ ] Add comprehensive docstrings to all `sc_tools/` modules
- [ ] Create API documentation

## 📝 Notes

- Library name: `sc-tools` (Python package: `sc_tools`)
- API follows scanpy pattern: `st.pl.*` for plotting, `st.tl.*` for tools
- All `sc_tools/` code should be generic and well-documented
- Scripts should import from `sc_tools.pl.*`, `sc_tools.tl.*`, etc.
- `scripts/dev/` scripts are temporary and should be cleaned up regularly
- `scripts/archive/` is read-only for reference only

## 🔄 Next Steps

1. (Optional) Refactor scripts to use `st.pl.spatial`, `st.pl.heatmaps`, `st.pl.save_figure` where applicable
2. Review and integrate imc-analysis functionalities
3. Move test scripts to `scripts/dev/`
4. Update existing scripts to use new API: `st.pl.*`, `st.tl.*`
5. Update `Architecture.md`
6. Project-specific modules (`sc_analysis/`) and script organization (phase folders, archive)
