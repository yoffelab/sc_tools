# ============================================================================
# GGO Visium Analysis Pipeline Makefile
# ============================================================================
# This Makefile orchestrates the complete spatial transcriptomics analysis
# pipeline from raw data ingestion to publication-ready figures.
#
# Usage:
#   make all              - Run entire pipeline
#   make phase1           - Run Phase I: Data Ingestion
#   make phase2           - Run Phase II: Preprocessing
#   make phase3           - Run Phase III: Gene Scoring & Deconvolution
#   make phase4           - Run Phase IV: Analysis & Statistics
#   make phase5           - Run Phase V: Visualization
#   make clean            - Remove intermediate results (use with caution)
#
# Dependencies: Python 3.x, conda environment (see environment.yml)
# ============================================================================

.PHONY: all phase1 phase2 phase3 phase4 phase5 clean help

# Default target: run entire pipeline
all: phase1 phase2 phase3 phase4 phase5
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
results/adata.annotation.h5ad: scripts/loupe2adata.py
	@echo "=========================================="
	@echo "Phase I.1: Converting Loupe data to AnnData"
	@echo "=========================================="
	python scripts/loupe2adata.py
	@echo "✅ Created: results/adata.annotation.h5ad"

# ----------------------------------------------------------------------------
# Script: annotation2mask_img.py
# Purpose: Add image masks and pathologist annotations to AnnData
# Input:  results/adata.annotation.h5ad
# Output: results/adata.annotation.masked.h5ad
# Notes:  Embeds tissue images and masks into AnnData object for spatial
#         visualization. Required for downstream spatial analysis.
# ----------------------------------------------------------------------------
results/adata.annotation.masked.h5ad: results/adata.annotation.h5ad scripts/annotation2mask_img.py
	@echo "=========================================="
	@echo "Phase I.2: Adding image masks and annotations"
	@echo "=========================================="
	python scripts/annotation2mask_img.py
	@echo "✅ Created: results/adata.annotation.masked.h5ad"

# Phase I target
phase1: results/adata.annotation.masked.h5ad
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
results/scvi.h5ad: results/adata.annotation.masked.h5ad scripts/preprocessing.py
	@echo "=========================================="
	@echo "Phase II.1: Preprocessing and batch correction"
	@echo "=========================================="
	python scripts/preprocessing.py
	@echo "✅ Created: results/scvi.h5ad"

# ----------------------------------------------------------------------------
# Script: clustering.py
# Purpose: Perform Leiden clustering on scVI latent space
# Input:  results/scvi.h5ad
# Output: results/scvi.leiden.h5ad
# Notes:  Clusters spots in the batch-corrected latent space for
#         downstream cell type annotation
# ----------------------------------------------------------------------------
results/scvi.leiden.h5ad: results/scvi.h5ad scripts/clustering.py
	@echo "=========================================="
	@echo "Phase II.2: Leiden clustering"
	@echo "=========================================="
	python scripts/clustering.py
	@echo "✅ Created: results/scvi.leiden.h5ad"

# ----------------------------------------------------------------------------
# Script: celltyping.py
# Purpose: Annotate clusters with cell type labels
# Input:  results/scvi.leiden.h5ad
# Output: results/scvi.leiden.phenotyped.h5ad
# Notes:  Adds cell type annotations to clusters based on marker genes
#         and reference signatures
# ----------------------------------------------------------------------------
results/scvi.leiden.phenotyped.h5ad: results/scvi.leiden.h5ad scripts/celltyping.py
	@echo "=========================================="
	@echo "Phase II.3: Cell type annotation"
	@echo "=========================================="
	python scripts/celltyping.py
	@echo "✅ Created: results/scvi.leiden.phenotyped.h5ad"

# Phase II target
phase2: results/scvi.leiden.phenotyped.h5ad
	@echo "✅ Phase II complete: Preprocessing and clustering"

# ============================================================================
# Phase III: Gene Signature Scoring & Cell Type Deconvolution
# ============================================================================
# Scores gene signatures and performs cell type deconvolution.

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
results/adata.img.genescores.h5ad: results/adata.annotation.masked.h5ad metadata/gene_signatures.json scripts/score_gene_signatures.py
	@echo "=========================================="
	@echo "Phase III.1: Gene signature scoring"
	@echo "=========================================="
	python scripts/score_gene_signatures.py
	@echo "✅ Created: results/adata.img.genescores.h5ad"

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
results/adata.deconvolution.h5ad: results/adata.img.genescores.h5ad results/seurat_object.h5ad scripts/celltype_deconvolution_phase3.py
	@echo "=========================================="
	@echo "Phase III.2: Cell type deconvolution"
	@echo "=========================================="
	python scripts/celltype_deconvolution_phase3.py
	@echo "✅ Created: results/adata.deconvolution.h5ad"

# Phase III target
phase3: results/adata.img.genescores.h5ad
	@echo "✅ Phase III complete: Gene scoring (deconvolution optional)"

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
results/tls_clustered.h5ad: results/adata.img.genescores.h5ad scripts/create_tls_anndata.py
	@echo "=========================================="
	@echo "Phase IV.1: TLS region extraction"
	@echo "=========================================="
	python scripts/create_tls_anndata.py
	@echo "✅ Created: results/tls_clustered.h5ad"

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
figures/manuscript/tls_analysis: results/adata.img.genescores.h5ad results/tls_clustered.h5ad scripts/tls_analysis.py
	@echo "=========================================="
	@echo "Phase IV.2: TLS analysis"
	@echo "=========================================="
	python scripts/tls_analysis.py
	@mkdir -p figures/manuscript/tls_analysis
	@echo "✅ Created: figures/manuscript/tls_analysis/"

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
figures/manuscript/spatial_analysis: results/adata.img.genescores.h5ad scripts/spatial_analysis.py
	@echo "=========================================="
	@echo "Phase IV.3: Spatial analysis (SVG and ligrec)"
	@echo "=========================================="
	python scripts/spatial_analysis.py
	@mkdir -p figures/manuscript/spatial_analysis
	@echo "✅ Created: figures/manuscript/spatial_analysis/"

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
figures/manuscript/tumor_differences: results/adata.img.genescores.h5ad scripts/tumor_differences.py
	@echo "=========================================="
	@echo "Phase IV.4: Tumor differences analysis"
	@echo "=========================================="
	python scripts/tumor_differences.py
	@mkdir -p figures/manuscript/tumor_differences
	@echo "✅ Created: figures/manuscript/tumor_differences/"

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
figures/manuscript/macrophage_localization: results/adata.img.genescores.h5ad scripts/macrophage_localization.py
	@echo "=========================================="
	@echo "Phase IV.5: Macrophage localization analysis"
	@echo "=========================================="
	python scripts/macrophage_localization.py
	@mkdir -p figures/manuscript/macrophage_localization
	@echo "✅ Created: figures/manuscript/macrophage_localization/"

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
figures/process_colocalization: results/adata.img.genescores.h5ad scripts/process_colocalization.py
	@echo "=========================================="
	@echo "Phase IV.6: Process colocalization analysis"
	@echo "=========================================="
	python scripts/process_colocalization.py
	@mkdir -p figures/process_colocalization
	@mkdir -p results/process_colocalization
	@echo "✅ Created: figures/process_colocalization/"

# Phase IV target
phase4: figures/manuscript/tumor_differences figures/manuscript/macrophage_localization figures/process_colocalization
	@echo "✅ Phase IV complete: Analysis and statistical testing"

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
figures/manuscript/spatial_plots: results/adata.img.genescores.h5ad scripts/manuscript_spatial_plots.py
	@echo "=========================================="
	@echo "Phase V.1: Publication-ready spatial plots"
	@echo "=========================================="
	python scripts/manuscript_spatial_plots.py
	@mkdir -p figures/manuscript/spatial_plots
	@echo "✅ Created: figures/manuscript/spatial_plots/"

# ----------------------------------------------------------------------------
# Script: signature_heatmap.py
# Purpose: Generate comprehensive heatmap and clustermap visualizations
# Input:  results/adata.img.genescores.h5ad
# Output: figures/manuscript/signature_heatmaps/ (4 figures total)
# Notes:  - Rows: Spots (sorted by patient and solidity)
#         - Columns: Gene signatures (clustered in clustermaps)
#         - Annotations: Patient and solidity color bars
#         - 3-level hierarchy: Macro grouping, secondary grouping, fine clustering
#         - 4 figures: 2 heatmaps (patient→solidity, solidity→patient) and 2 clustermaps
#         - Solidity legend in top right corner
# ----------------------------------------------------------------------------
figures/manuscript/signature_heatmaps: results/adata.img.genescores.h5ad scripts/signature_heatmap.py
	@echo "=========================================="
	@echo "Phase V.2: Signature heatmap visualization"
	@echo "=========================================="
	python scripts/signature_heatmap.py
	@mkdir -p figures/manuscript/signature_heatmaps
	@echo "✅ Created: figures/manuscript/signature_heatmaps/"

# Phase V target
phase5: figures/manuscript/spatial_plots figures/manuscript/signature_heatmaps
	@echo "✅ Phase V complete: Visualization"

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
	@echo "  make phase3           - Phase III: Gene scoring and deconvolution"
	@echo "  make phase4           - Phase IV: Analysis and statistics"
	@echo "  make phase5           - Phase V: Visualization"
	@echo "  make clean            - Remove intermediate results (CAUTION)"
	@echo "  make help             - Show this help message"
	@echo ""
	@echo "Individual script targets:"
	@echo "  results/adata.annotation.masked.h5ad     - Data ingestion"
	@echo "  results/scvi.leiden.phenotyped.h5ad      - Preprocessing"
	@echo "  results/adata.img.genescores.h5ad        - Gene scoring"
	@echo "  results/adata.deconvolution.h5ad          - Deconvolution"
	@echo "  results/tls_clustered.h5ad                - TLS extraction"
	@echo "  figures/manuscript/tumor_differences      - Tumor analysis"
	@echo "  figures/manuscript/macrophage_localization - Macrophage analysis"
	@echo "  figures/process_colocalization            - Process colocalization"
	@echo "  figures/manuscript/signature_heatmaps      - Signature heatmaps"

# Clean target (use with caution - removes intermediate files)
clean:
	@echo "⚠️  WARNING: This will remove intermediate results"
	@echo "Press Ctrl+C to cancel, or Enter to continue..."
	@read
	rm -f results/*.h5ad
	rm -f results/raw.concat.h5ad
	rm -rf figures/manuscript/*
	rm -rf figures/stats/*
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

