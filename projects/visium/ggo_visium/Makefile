# ============================================================================
# Spatial Omics Analysis Pipeline Makefile
# ============================================================================
# Project-aware: set PROJECT to the project dir (e.g. projects/visium/ggo_visium).
# Scripts run with cwd = PROJECT so paths (results/, metadata/, etc.) are correct.
#
# Usage:
#   make all              - Run entire pipeline (default PROJECT=projects/visium/ggo_visium)
#   make PROJECT=projects/visium/ggo_visium phase1
#   make phase1           - Run Phase I: Data Ingestion
#   make phase2           - Run Phase II: Preprocessing
#   ... phase3, phase4, phase5
#   make clean            - Remove intermediate results (use with caution)
#
# Dependencies: Python 3.x, conda environment. Activate before running:
#   conda activate robin
# ============================================================================

PROJECT ?= projects/visium/ggo_visium
CONDA_ENV ?= robin

.PHONY: all phase1 phase2 phase3 phase3.5b phase4 phase5 phase4-direct phase5-direct signature_heatmaps_versioned clean help

# Default target: run entire pipeline
all: phase1 phase2 phase3.5b phase4 phase5
	@echo "=========================================="
	@echo "✅ Complete pipeline finished successfully"
	@echo "=========================================="

# ============================================================================
# Phase I: Data Ingestion & Annotation
# ============================================================================
# Converts raw spatial transcriptomics data to standardized AnnData format
# and adds pathologist annotations and image masks.

# ----------------------------------------------------------------------------
# Script: loupe2adata.py
# Purpose: Convert 10x Genomics Loupe Browser exports to AnnData format
# Input:  Raw spatial data files (from data/ directory)
# Output: results/adata.annotation.h5ad
#         results/raw.concat.h5ad
# Notes:  First step in pipeline - creates base AnnData object with
#         spatial coordinates and basic metadata
# ----------------------------------------------------------------------------
$(PROJECT)/results/adata.annotation.h5ad: $(PROJECT)/scripts/loupe2adata.py
	@echo "=========================================="
	@echo "Phase I.1: Converting Loupe data to AnnData"
	@echo "=========================================="
	cd $(PROJECT) && python scripts/loupe2adata.py
	@echo "✅ Created: $(PROJECT)/results/adata.annotation.h5ad"

# ----------------------------------------------------------------------------
# Script: annotation2mask_img.py
# Purpose: Add image masks and pathologist annotations to AnnData
# Input:  results/adata.annotation.h5ad
# Output: results/adata.annotation.masked.h5ad
# Notes:  Embeds tissue images and masks into AnnData object for spatial
#         visualization. Required for downstream spatial analysis.
# ----------------------------------------------------------------------------
$(PROJECT)/results/adata.annotation.masked.h5ad: $(PROJECT)/results/adata.annotation.h5ad $(PROJECT)/scripts/annotation2mask_img.py
	@echo "=========================================="
	@echo "Phase I.2: Adding image masks and annotations"
	@echo "=========================================="
	cd $(PROJECT) && python scripts/annotation2mask_img.py
	@echo "✅ Created: $(PROJECT)/results/adata.annotation.masked.h5ad"

# Phase I target
phase1: $(PROJECT)/results/adata.annotation.masked.h5ad
	@echo "✅ Phase I complete: Data ingestion and annotation"

# ============================================================================
# Phase II: Preprocessing & Integration
# ============================================================================
# Performs quality control, normalization, batch correction, and clustering.

# ----------------------------------------------------------------------------
# Script: preprocessing.py
# Purpose: QC, normalization, and scVI batch correction
# Input:  results/adata.annotation.masked.h5ad
# Output: results/scvi.h5ad
# Notes:  - Filters cells and genes based on QC metrics
#         - Normalizes and log-transforms data
#         - Trains scVI model for batch correction and latent representation
#         - Calculates mitochondrial, ribosomal, and hemoglobin content
#         - Identifies highly variable genes (HVG) and spatially variable genes (SVG)
# ----------------------------------------------------------------------------
$(PROJECT)/results/scvi.h5ad: $(PROJECT)/results/adata.annotation.masked.h5ad $(PROJECT)/scripts/preprocessing.py
	@echo "=========================================="
	@echo "Phase II.1: Preprocessing and batch correction"
	@echo "=========================================="
	cd $(PROJECT) && python scripts/preprocessing.py
	@echo "✅ Created: $(PROJECT)/results/scvi.h5ad"

# ----------------------------------------------------------------------------
# Script: clustering.py
# Purpose: Perform Leiden clustering on scVI latent space
# Input:  results/scvi.h5ad
# Output: results/scvi.leiden.h5ad
# Notes:  Clusters spots in the batch-corrected latent space.
#         Automated cell typing is in Phase 3.5b.
# ----------------------------------------------------------------------------
$(PROJECT)/results/scvi.leiden.h5ad: $(PROJECT)/results/scvi.h5ad $(PROJECT)/scripts/clustering.py
	@echo "=========================================="
	@echo "Phase II.2: Leiden clustering"
	@echo "=========================================="
	cd $(PROJECT) && python scripts/clustering.py
	@echo "✅ Created: $(PROJECT)/results/scvi.leiden.h5ad"

# Phase II target (preprocessing + clustering only; no cell typing)
phase2: $(PROJECT)/results/scvi.leiden.h5ad
	@echo "✅ Phase II complete: Preprocessing and clustering"

# ============================================================================
# Phase 3.5b: Gene Signature Scoring, Automated Cell Typing, Deconvolution
# ============================================================================
# Separate branch from Phase 3; connects to Phase 4 (Phase 4 skippable if adequate).

# ----------------------------------------------------------------------------
# Script: celltyping.py
# Purpose: Annotate clusters with cell type labels (automated)
# Input:  results/scvi.leiden.h5ad
# Output: results/scvi.leiden.phenotyped.h5ad
# Notes:  Adds cell type annotations to clusters based on marker genes
#         and reference signatures
# ----------------------------------------------------------------------------
$(PROJECT)/results/scvi.leiden.phenotyped.h5ad: $(PROJECT)/results/scvi.leiden.h5ad $(PROJECT)/scripts/celltyping.py
	@echo "=========================================="
	@echo "Phase 3.5b.0: Automated cell typing"
	@echo "=========================================="
	cd $(PROJECT) && python scripts/celltyping.py
	@echo "✅ Created: $(PROJECT)/results/scvi.leiden.phenotyped.h5ad"

# ----------------------------------------------------------------------------
# Script: score_gene_signatures.py
# Purpose: Calculate and store gene signature scores in AnnData
# Input:  results/adata.annotation.masked.h5ad
#         metadata/gene_signatures.json
# Output: results/adata.img.genescores.h5ad
# Notes:  - Scores all gene signatures from metadata/gene_signatures.json
#         - Includes tumor programs (EMT, Hypoxia, Proliferative)
#         - Includes immune programs (Lymphoid, Innate, Processes)
#         - Includes Liron signatures (T cells, macrophages)
#         - Z-scores each signature across spots for visualization
#         - Required for all downstream analysis scripts
# ----------------------------------------------------------------------------
$(PROJECT)/results/adata.img.genescores.h5ad: $(PROJECT)/results/adata.annotation.masked.h5ad $(PROJECT)/metadata/gene_signatures.json $(PROJECT)/scripts/score_gene_signatures.py
	@echo "=========================================="
	@echo "Phase 3.5b.1: Gene signature scoring"
	@echo "=========================================="
	cd $(PROJECT) && python scripts/score_gene_signatures.py
	@echo "✅ Created: $(PROJECT)/results/adata.img.genescores.h5ad"

# ----------------------------------------------------------------------------
# Script: celltype_deconvolution_phase3.py
# Purpose: Deconvolve cell type proportions in spatial spots
# Input:  results/seurat_object.h5ad (single-cell reference)
#         results/adata.img.genescores.h5ad (spatial data)
# Output: results/adata.deconvolution.h5ad
# Notes:  - Uses multiple methods with fallback: DestVI, Cell2location, Tangram
#         - Processes libraries individually to manage memory
#         - Includes comprehensive memory profiling and optimization
#         - Estimates cell type proportions for each spatial spot
#         - Required for differential cell-type proportion analysis
# ----------------------------------------------------------------------------
$(PROJECT)/results/adata.deconvolution.h5ad: $(PROJECT)/results/adata.img.genescores.h5ad $(PROJECT)/results/seurat_object.h5ad $(PROJECT)/scripts/celltype_deconvolution_phase3.py
	@echo "=========================================="
	@echo "Phase 3.5b.2: Cell type deconvolution"
	@echo "=========================================="
	cd $(PROJECT) && python scripts/celltype_deconvolution_phase3.py
	@echo "✅ Created: $(PROJECT)/results/adata.deconvolution.h5ad"

# Phase 3.5b target (automated cell typing + gene scoring required; deconvolution optional)
phase3.5b: $(PROJECT)/results/scvi.leiden.phenotyped.h5ad $(PROJECT)/results/adata.img.genescores.h5ad
	@echo "✅ Phase 3.5b complete: Cell typing + gene scoring (deconvolution optional)"

# Backward compatibility: phase3 runs phase 3.5b
phase3: phase3.5b

# ============================================================================
# Phase IV: Analysis & Statistical Testing
# ============================================================================
# Performs statistical comparisons and specialized analyses.

# ----------------------------------------------------------------------------
# Script: create_tls_anndata.py
# Purpose: Extract and aggregate TLS regions from all tissue slides
# Input:  results/adata.img.genescores.h5ad
# Output: results/tls_clustered.h5ad
# Notes:  - Identifies TLS regions based on architecture annotations
#         - Creates aggregated AnnData object for TLS-specific analysis
#         - Required for tls_analysis.py
# ----------------------------------------------------------------------------
$(PROJECT)/results/tls_clustered.h5ad: $(PROJECT)/results/adata.img.genescores.h5ad $(PROJECT)/scripts/create_tls_anndata.py
	@echo "=========================================="
	@echo "Phase IV.1: TLS region extraction"
	@echo "=========================================="
	cd $(PROJECT) && python scripts/create_tls_anndata.py
	@echo "✅ Created: $(PROJECT)/results/tls_clustered.h5ad"

# ----------------------------------------------------------------------------
# Script: tls_analysis.py
# Purpose: TLS-specific transcriptomics analysis
# Input:  results/adata.img.genescores.h5ad
#         results/tls_clustered.h5ad
# Output: figures/manuscript/tls_analysis/ (various plots)
#         figures/stats/tls_*.csv (statistics)
# Notes:  - UMAP and spatial plots of TLSs
#         - B-cell/T-cell state comparison
#         - Ligand-receptor crosstalk analysis
#         - Statistical testing with FDR correction
# ----------------------------------------------------------------------------
$(PROJECT)/figures/manuscript/tls_analysis: $(PROJECT)/results/adata.img.genescores.h5ad $(PROJECT)/results/tls_clustered.h5ad $(PROJECT)/scripts/tls_analysis.py
	@echo "=========================================="
	@echo "Phase IV.2: TLS analysis"
	@echo "=========================================="
	cd $(PROJECT) && python scripts/tls_analysis.py
	@mkdir -p $(PROJECT)/figures/manuscript/tls_analysis
	@echo "✅ Created: $(PROJECT)/figures/manuscript/tls_analysis/"

# ----------------------------------------------------------------------------
# Script: spatial_analysis.py
# Purpose: Spatially Variable Genes (SVG) and ligand-receptor analysis
# Input:  results/adata.img.genescores.h5ad
# Output: figures/manuscript/spatial_analysis/ (SVG plots)
#         figures/stats/spatial_*.csv (SVG statistics)
# Notes:  - Identifies SVGs using Moran's I or squidpy.gr.spatial_autocorr
#         - Compares SVG expression between Non-Solid and Solid tumors
#         - Infers ligand-receptor interactions using squidpy.gr.ligrec
# ----------------------------------------------------------------------------
$(PROJECT)/figures/manuscript/spatial_analysis: $(PROJECT)/results/adata.img.genescores.h5ad $(PROJECT)/scripts/spatial_analysis.py
	@echo "=========================================="
	@echo "Phase IV.3: Spatial analysis (SVG and ligrec)"
	@echo "=========================================="
	cd $(PROJECT) && python scripts/spatial_analysis.py
	@mkdir -p $(PROJECT)/figures/manuscript/spatial_analysis
	@echo "✅ Created: $(PROJECT)/figures/manuscript/spatial_analysis/"

# ----------------------------------------------------------------------------
# Script: tumor_differences.py
# Purpose: 1-vs-all statistical testing across tumor types
# Input:  results/adata.img.genescores.h5ad
# Output: figures/manuscript/tumor_differences/ (boxplots)
#         figures/stats/tumor_differences_1vsall_stats.csv
# Notes:  - Tests tumor programs (EMT, Hypoxia, Proliferative)
#         - Tests immune programs (Immune_Lymphoid, Innate_Other, Processes, Liron)
#         - Performs Mann-Whitney U test with Benjamini-Hochberg correction
#         - Generates boxplots with significance bars and p-values
#         - One line per comparison in text annotations
# ----------------------------------------------------------------------------
$(PROJECT)/figures/manuscript/tumor_differences: $(PROJECT)/results/adata.img.genescores.h5ad $(PROJECT)/scripts/tumor_differences.py
	@echo "=========================================="
	@echo "Phase IV.4: Tumor differences analysis"
	@echo "=========================================="
	cd $(PROJECT) && python scripts/tumor_differences.py
	@mkdir -p $(PROJECT)/figures/manuscript/tumor_differences
	@echo "✅ Created: $(PROJECT)/figures/manuscript/tumor_differences/"

# ----------------------------------------------------------------------------
# Script: macrophage_localization.py
# Purpose: Measure macrophage colocalization with proliferative program
# Input:  results/adata.img.genescores.h5ad
# Output: figures/manuscript/macrophage_localization/ (scatterplots, heatmaps)
#         figures/stats/macrophage_proliferation_correlations.csv
# Notes:  - Stratifies by Normal, Non-Solid, and Solid tumors
#         - X-axis: Proliferative score, Y-axis: Macrophage scores
#         - Filters spots with proliferative < 0 or macrophage < 0 (gray)
#         - Calculates Pearson correlation per tumor type
#         - Includes regression lines for each condition
#         - Annotates correlations outside plot
#         - Generates correlation heatmap and clustermap
# ----------------------------------------------------------------------------
$(PROJECT)/figures/manuscript/macrophage_localization: $(PROJECT)/results/adata.img.genescores.h5ad $(PROJECT)/scripts/macrophage_localization.py
	@echo "=========================================="
	@echo "Phase IV.5: Macrophage localization analysis"
	@echo "=========================================="
	cd $(PROJECT) && python scripts/macrophage_localization.py
	@mkdir -p $(PROJECT)/figures/manuscript/macrophage_localization
	@echo "✅ Created: $(PROJECT)/figures/manuscript/macrophage_localization/"

# ----------------------------------------------------------------------------
# Script: process_colocalization.py
# Purpose: Analyze spatial co-occurrence patterns of processes across spots
# Input:  results/adata.img.genescores.h5ad
# Output: figures/process_colocalization/ (heatmaps, volcano plots)
#         results/process_colocalization/*.csv (statistics)
# Notes:  - Pearson correlation of processes across spots
#         - Moran's I for spatial autocorrelation
#         - Neighborhood enrichment (optional, disabled by default)
#         - Faceted volcano plots for Normal vs Non-Solid vs Solid comparisons
#         - Signature filtering (include/exclude) for selective visualization
#         - Correlation heatmap and clustermap
# ----------------------------------------------------------------------------
$(PROJECT)/figures/process_colocalization: $(PROJECT)/results/adata.img.genescores.h5ad $(PROJECT)/scripts/process_colocalization.py
	@echo "=========================================="
	@echo "Phase IV.6: Process colocalization analysis"
	@echo "=========================================="
	cd $(PROJECT) && python scripts/process_colocalization.py
	@mkdir -p $(PROJECT)/figures/process_colocalization
	@mkdir -p $(PROJECT)/results/process_colocalization
	@echo "✅ Created: $(PROJECT)/figures/process_colocalization/"

# Phase IV target
phase4: $(PROJECT)/figures/manuscript/tumor_differences $(PROJECT)/figures/manuscript/macrophage_localization $(PROJECT)/figures/process_colocalization
	@echo "✅ Phase IV complete: Analysis and statistical testing"

# Phase IV direct: run phase 4 scripts only (assumes adata.img.genescores.h5ad exists).
# Use when phase 1-3 results exist and you want to skip dependency checks.
phase4-direct:
	@echo "=========================================="
	@echo "Phase IV (direct): Running analysis scripts"
	@echo "=========================================="
	cd $(PROJECT) && python scripts/tumor_differences.py
	cd $(PROJECT) && python scripts/macrophage_localization.py
	cd $(PROJECT) && python scripts/process_colocalization.py
	@mkdir -p $(PROJECT)/figures/manuscript/tumor_differences $(PROJECT)/figures/manuscript/macrophage_localization $(PROJECT)/figures/process_colocalization $(PROJECT)/results/process_colocalization
	@echo "✅ Phase IV (direct) complete"

# ============================================================================
# Phase V: Visualization & Publication Figures
# ============================================================================
# Generates publication-ready figures following statistical annotation standards.

# ----------------------------------------------------------------------------
# Script: manuscript_spatial_plots.py
# Purpose: Generate publication-ready spatial plots
# Input:  results/adata.img.genescores.h5ad
# Output: figures/manuscript/ (various spatial visualizations)
# Notes:  - Spatial plots of gene signatures and cell types
#         - Follows statistical annotation rules from skills.md
#         - All plots include significance bars and p-values
# ----------------------------------------------------------------------------
$(PROJECT)/figures/manuscript/spatial_plots: $(PROJECT)/results/adata.img.genescores.h5ad $(PROJECT)/scripts/manuscript_spatial_plots.py
	@echo "=========================================="
	@echo "Phase V.1: Publication-ready spatial plots"
	@echo "=========================================="
	cd $(PROJECT) && python scripts/manuscript_spatial_plots.py
	@mkdir -p $(PROJECT)/figures/manuscript/spatial_plots
	@echo "✅ Created: $(PROJECT)/figures/manuscript/spatial_plots/"

# ----------------------------------------------------------------------------
# Script: signature_heatmap.py (legacy)
# Purpose: Generate heatmap/clustermap visualizations (no datetime stamp)
# ----------------------------------------------------------------------------
$(PROJECT)/figures/manuscript/signature_heatmaps: $(PROJECT)/results/adata.img.genescores.h5ad $(PROJECT)/scripts/signature_heatmap.py
	@echo "=========================================="
	@echo "Phase V.2: Signature heatmap visualization (legacy)"
	@echo "=========================================="
	cd $(PROJECT) && python scripts/signature_heatmap.py
	@mkdir -p $(PROJECT)/figures/manuscript/signature_heatmaps
	@echo "✅ Created: $(PROJECT)/figures/manuscript/signature_heatmaps/"

# ----------------------------------------------------------------------------
# Script: signature_heatmap_versioned.py (preferred: sc_tools + datetime stamp)
# Purpose: Same heatmaps/clustermaps using sc_tools.pl; outputs versioned PDF/PNG
# Output: figures/manuscript/signature_heatmaps/pdf/YYDDMM.hh.mm.<basename>.pdf
#         figures/manuscript/signature_heatmaps/png/YYDDMM.hh.mm.<basename>.png
# ----------------------------------------------------------------------------
signature_heatmaps_versioned:
	@echo "=========================================="
	@echo "Signature heatmaps (versioned, sc_tools)"
	@echo "=========================================="
	cd $(PROJECT) && python scripts/signature_heatmap_versioned.py
	@echo "✅ Versioned heatmaps in $(PROJECT)/figures/manuscript/signature_heatmaps/pdf and .../png"

# Phase V target
phase5: $(PROJECT)/figures/manuscript/spatial_plots $(PROJECT)/figures/manuscript/signature_heatmaps
	@echo "✅ Phase V complete: Visualization"

# Phase V direct: run phase 5 scripts only (assumes adata.img.genescores.h5ad exists).
# Uses versioned signature heatmap script (datetime-stamped figures).
phase5-direct:
	@echo "=========================================="
	@echo "Phase V (direct): Running visualization scripts"
	@echo "=========================================="
	cd $(PROJECT) && python scripts/manuscript_spatial_plots.py
	cd $(PROJECT) && python scripts/signature_heatmap_versioned.py
	@mkdir -p $(PROJECT)/figures/manuscript/spatial_plots $(PROJECT)/figures/manuscript/signature_heatmaps
	@echo "✅ Phase V (direct) complete"

# ============================================================================
# Utility Targets
# ============================================================================

# Help target
help:
	@echo "GGO Visium Analysis Pipeline Makefile"
	@echo ""
	@echo "Available targets:"
	@echo "  make all              - Run entire pipeline (Phases I-V)"
	@echo "  make phase1           - Phase I: Data ingestion and annotation"
	@echo "  make phase2           - Phase II: Preprocessing and clustering"
	@echo "  make phase3.5b        - Phase 3.5b: Gene scoring and deconvolution"
	@echo "  make phase3           - Alias for phase3.5b"
	@echo "  make phase4           - Phase IV: Analysis and statistics"
	@echo "  make phase5           - Phase V: Visualization"
	@echo "  make phase4-direct    - Phase IV only (assumes genescores exist; use conda activate $(CONDA_ENV))"
	@echo "  make phase5-direct    - Phase V only (versioned heatmap output)"
	@echo "  make signature_heatmaps_versioned - Signature heatmaps only (datetime-stamped pdf/png)"
	@echo "  make clean            - Remove intermediate results (CAUTION)"
	@echo "  make help             - Show this help message"
	@echo ""
	@echo "Current PROJECT: $(PROJECT)"
	@echo "Individual script targets (paths under PROJECT):"
	@echo "  results/adata.annotation.masked.h5ad       - Data ingestion"
	@echo "  results/scvi.leiden.phenotyped.h5ad        - Phase 3.5b (automated cell typing)"
	@echo "  results/adata.img.genescores.h5ad           - Gene scoring (Phase 3.5b)"
	@echo "  results/adata.deconvolution.h5ad          - Deconvolution"
	@echo "  results/tls_clustered.h5ad                 - TLS extraction"
	@echo "  figures/manuscript/tumor_differences       - Tumor analysis"
	@echo "  figures/manuscript/macrophage_localization - Macrophage analysis"
	@echo "  figures/process_colocalization             - Process colocalization"
	@echo "  figures/manuscript/signature_heatmaps     - Signature heatmaps (versioned: pdf/ and png/ with YYDDMM.hh.mm prefix)"

# Clean target (use with caution - removes intermediate files)
clean:
	@echo "⚠️  WARNING: This will remove intermediate results for $(PROJECT)"
	@echo "Press Ctrl+C to cancel, or Enter to continue..."
	@read
	rm -f $(PROJECT)/results/*.h5ad
	rm -f $(PROJECT)/results/raw.concat.h5ad
	rm -rf $(PROJECT)/figures/manuscript/*
	rm -rf $(PROJECT)/figures/stats/*
	@echo "✅ Cleaned intermediate results"

# ============================================================================
# Notes for Maintenance
# ============================================================================
# When adding new scripts:
# 1. Add a new target with appropriate dependencies
# 2. Document the script's purpose, inputs, and outputs
# 3. Add to appropriate phase target
# 4. Update help target if needed
#
# Script execution order is critical - ensure dependencies are correct!
# All scripts should follow the standards in skills.md for statistical rigor.
# ============================================================================

