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
  - ✅ `qc/` - QC module (placeholder; planned: metrics, spatial, plots per WORKFLOW.md)
- ✅ `projects/` - Scalable project layout (replaces top-level sc_analysis and root scripts/metadata)
  - ✅ `projects/visium/`, `visium_hd/`, `xenium/`, `imc/` - Data-type folders
  - ✅ `projects/visium/ggo_visium/` - Example project with data/, figures/, metadata/, scripts/, results/, outputs/
  - ✅ `create_project.sh` - Creates new project dirs; usage: `./create_project.sh <project_name> <data_type>`
  - ✅ `migrate_to_ggo_visium.sh` - One-time move of root scripts/metadata into projects/visium/ggo_visium/
- ✅ Makefile project-aware: `PROJECT ?= projects/visium/ggo_visium`; all paths use `$(PROJECT)/...`

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

## 🚧 In Progress — Next immediate steps

**0. Testing (implementation order: 1st ggo_visium, 2nd sc_tools, 3rd functions)**
- [ ] **ggo_visium project tests:** Create `projects/visium/ggo_visium/tests/`. Integration/smoke tests for Makefile and scripts. Phases 1–3 not fully tested; Phases 4–5 should work (or with few fixes).
- [ ] **sc_tools package tests:** Create `sc_tools/tests/`. Unit tests for pl, tl, qc with synthetic fixtures.
- [ ] All new code must compile and pass tests before merge.

**1. Align Makefile with new 7-phase workflow**
- [ ] Update Makefile targets to match Phases 1–7 (WORKFLOW.md).
- [ ] Verify Phase 3–5 (preprocessing, downstream) run successfully.

**1b. Script path migration (metadata move)**
- [ ] Update scripts to use $(PROJECT)/metadata/ instead of metadata/. Affected: score_gene_signatures.py, tumor_differences.py, tls_analysis.py, manuscript_spatial_plots.py, celltyping.py, etc.
- [ ] Confirm outputs match previous results (versioned filenames acceptable).
- [ ] Fix any broken targets or paths.

**2. Modular scripts: config + sc_tools imports**
- [ ] Refactor scripts to use config files (dicts of arguments) and import from `sc_tools` directly.
- [ ] No duplicated logic in scripts; thin orchestration only (config → st.pl / st.tl / st.data → versioned outputs).

## 📋 To Do (later)

### 5. Project-Specific Modules
- Project-specific code lives under `projects/<data_type>/<project_name>/scripts/`.
- [ ] Optional: Extract shared logic (TLS, macrophage, signatures) into `sc_tools` when reused across data types.

### 6. Script Organization
- [x] Scalable layout: `projects/visium|visium_hd|xenium|imc/<name>/` with data, figures, metadata, scripts, results, outputs.
- [x] Single legacy folder per project: `scripts/old_code/` inside each project.
- [ ] Organize production scripts by phase within project if desired.
- [ ] Update imports in all scripts to use `sc_tools.pl.*`, `sc_tools.tl.*`

### 7. Integration with imc-analysis
- [ ] Review imc-analysis library structure
- [ ] Identify relevant functionalities
- [ ] Adapt and integrate into `sc_tools`

### 8. Testing
- [ ] Add pytest to requirements / pyproject.toml
- [ ] Document test run commands in README / Architecture

### 9. Documentation
- [ ] Update `Architecture.md` with new structure
- [ ] Create migration guide for updating existing scripts
- [ ] Add comprehensive docstrings to all `sc_tools/` modules
- [ ] Create API documentation

## 📝 Notes

- Library name: `sc-tools` (Python package: `sc_tools`)
- API follows scanpy pattern: `st.pl.*` for plotting, `st.tl.*` for tools
- All `sc_tools/` code should be generic and well-documented
- Scripts should import from `sc_tools.pl.*`, `sc_tools.tl.*`, etc.
- New projects: `./create_project.sh <project_name> visium|visium_hd|xenium|imc`. Run make with `PROJECT=projects/<type>/<name>` or use default `projects/visium/ggo_visium`.
- Legacy: `scripts/old_code/` inside each project is read-only for reference.

## 🔄 Next Steps (after immediate steps)

1. Move test scripts to `scripts/dev/`, archive to `scripts/archive/`
2. Project-specific modules (`sc_analysis/`) and script organization by phase
3. Review and integrate imc-analysis functionalities
4. Update `Architecture.md` and add API documentation
