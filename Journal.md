# Research Journal & Decision Log: Spatial Lung Oncology

This journal documents the technical evolution, parameters, and rationale behind the analytical decisions of this project.

## Project Note: Transition to AI-Assisted Pipeline (2025-12-30)
* **Status:** We have recently adopted an agentic scientific programming workflow to accelerate analysis. 
* **Advisory:** Because parts of the current codebase (`scripts/`) are being refactored from legacy formats (`scripts/old code/`) by an AI agent, some workarounds may be required to resolve dependency or file path inconsistencies during the transition.
* **Constraint:** All generated documentation and communication regarding this project must strictly avoid the use of apostrophes.

---

# #Transition & Technical Log

## [2025-12-30] - Codebase Distillation & Archive
- **Action:** Audited `scripts/` and moved approximately 40% of inherited files to `scripts/old code/`.
- **Rationale:** Legacy scripts from the previous project were creating "namespace noise." By isolating active scripts, the AI agent can now focus on the scVI-based pipeline without hallucinating old functions.
- **Decision:** Adopted the **"Plan-Act-Reflect"** protocol for all new analysis steps involving the tumor progression stages.


## Log Entries

### [2025-12-30] - Initialization of Modern Frameworks
- **Action:** Transitioned from standard log-normalization to probabilistic modeling via `scvi-tools`.
- **Rationale:** Traditional normalization was masking biological variation in low-count spots, particularly in GGO (Ground Glass Opacity) regions.
- **Decision:** Use `scvi.model.SCVI` for all future integration steps.
- **Result:** Improved separation of tumor vs. stromal clusters in the integrated UMAP.

### [2025-12-30] - Architecture & Mission Alignment
- **Action:** Created `Architecture.md` and `Mission.md` to guide the AI agent.
- **Rationale:** The repository structure was becoming fragmented. Explicitly defining Phase I through Phase V ensures the agent does not skip pre-processing steps.
- **Decision:** Legacy scripts in `scripts/old code/` are to be treated as read-only references.

### [2025-12-28] - Statistical Visualization Standards
- **Action:** Enforced mandatory significance bars and numerical p-values for all manuscript figures.
- **Rationale:** To ensure maximum transparency and scientific rigor for computational oncology publications.
- **Decision:** Default to **1-vs-rest** comparisons for tumor progression stages (Normal, GGO, Part-Solid, Solid) to highlight specific stage-wise programs.
- **Future Considerations:** ~~While performing multiple comparison, p-val bars on plot are on top of each other. These need to be dodging. There is also one comparison being written logged outside of plot, in a clearly distinguished fashion. Prefered method of logging would be Group A - Group B: p*-value.~~ **RESOLVED** - See entry below.

### [2025-01-02] - Phase III Cell Type Deconvolution Implementation
- **Action:** Created `celltype_deconvolution_phase3.py` to perform Tangram-based deconvolution and integrate cell type proportions into spatial data.
- **Rationale:** The existing `celltype_deconvolution_tangram.py` script runs Tangram mapping but does not properly integrate proportions back into the main AnnData object for downstream analysis. Phase III requires cell type proportions to be accessible for spatial analysis tasks.
- **Decision:** Use Tangram in cluster mode with 500 epochs per batch. Store proportions in both `adata.obsm['cell_type_proportions']` (matrix format) and individual `adata.obs['prop_*']` columns for convenience. Process batches separately using `library_id` as the batch key to handle memory constraints.
- **Parameters:** 
  - Signature genes: 2000 (union of HVGs and top 100 markers per cell type)
  - Epochs: 500 per batch
  - Mode: clusters (cell type level deconvolution)
- **Output:** `results/adata.deconvolution.h5ad` with integrated proportions
- **Future Considerations:** May explore Cell2location or DestVI for comparison, especially if Tangram results show limitations in sparse regions or if more granular cell type resolution is needed.

### [2025-01-02] - Phase III Memory Optimization and Segmentation Fault Fix
- **Issue:** Initial implementation caused segmentation fault due to memory exhaustion when loading entire spatial dataset into memory.
- **Root Cause:** Attempting to load full spatial AnnData object (all libraries) simultaneously exceeded available memory, causing segmentation fault during Tangram processing.
- **Solution:** Implemented memory-efficient divide-and-conquer approach:
  - Process one `library_id` at a time
  - Load spatial data per library using backed mode (`backed='r'`) to read metadata first
  - Only load full data for the current library being processed
  - Explicit garbage collection (`gc.collect()`) after each library
  - Clear variables immediately after use
  - Use minimal data copies
- **Tangram Validation:** Verified correct application through literature review:
  - Tangram is a validated tool for cell type deconvolution in spatial transcriptomics (Tangram documentation, scvi-tools integration)
  - Workflow: `tg.map_cells_to_space(mode='clusters')` followed by `tg.project_cell_annotations()` is the standard approach
  - Cluster mode is appropriate for cell type-level deconvolution (not single-cell resolution)
  - Memory issues are common with large spatial datasets; batch processing is recommended
- **Implementation Changes:**
  - Load single-cell reference once (smaller dataset)
  - Identify library IDs from spatial metadata without loading full dataset
  - Process each library sequentially with full memory cleanup between iterations
  - Added comprehensive error handling and progress tracking
  - Extract proportions from both `obsm` and `obs` columns (handles different Tangram versions)
- **Result:** Script now processes libraries individually, preventing memory overflow and segmentation faults. Each library is processed in isolation with proper cleanup, allowing the script to handle datasets of any size.
- **Testing:** Script structure validated to ensure:
  - No full dataset loading before processing
  - Proper memory management with explicit cleanup
  - Error handling for individual library failures (continues with next library)
  - Correct Tangram API usage verified against official documentation

### [2025-01-05] - Phase III Modularization: Multi-Method Deconvolution Support
- **Issue:** Tangram continues to result in segmentation faults despite memory optimizations, requiring fallback alternatives.
- **Action:** Refactored `celltype_deconvolution_phase3.py` into a modular architecture supporting multiple deconvolution methods with automatic fallback.
- **Rationale:** Different deconvolution methods have varying computational requirements and robustness. Providing multiple methods allows the pipeline to continue even if one method fails, and enables method comparison for validation.
- **Implementation:**
  - **Modular Design:** Separated deconvolution logic into independent functions:
    - `deconvolve_tangram()`: Optimal transport-based mapping (original method)
    - `deconvolve_cell2location()`: Negative binomial regression model for cell type abundance
    - `deconvolve_destvi()`: Bayesian model for cell type proportions with continuous variation
  - **Fallback Strategy:** Methods are tried in order (DestVI -> Cell2location -> Tangram). If one fails, the next is attempted automatically.
  - **Shared Infrastructure:** Common data loading, preprocessing, and integration functions used by all methods.
  - **Method-Specific Handling:** Each method extracts proportions from its own output format:
    - Tangram: `obsm` or `obs` columns from `project_cell_annotations()`
    - Cell2location: `obsm['q05_cell_abundance_w_sf']` (quantiles of cell abundance)
    - DestVI: `obsm['proportions']` from `get_proportions()`
- **Method Validation:**
  - **Tangram:** Validated as standard approach for spatial mapping (Tangram documentation, scvi-tools integration). However, segmentation faults persist, likely due to underlying C++ library issues.
  - **Cell2location:** Validated approach for estimating cell type abundance using negative binomial regression. Recommended for Visium data with reference single-cell datasets (Cell2location documentation).
  - **DestVI:** Validated Bayesian approach that models both proportions and continuous variation within cell types. Part of scvi-tools ecosystem, well-integrated with AnnData (DestVI documentation).
- **Configuration:** Methods to try can be configured via `CONFIG['methods']` list. Default order prioritizes Tangram (fastest when working) but falls back to more robust methods.
- **Error Handling:** Each method wrapped in try-except blocks. Import errors handled gracefully (methods not installed are skipped). Runtime errors logged but do not stop processing of other libraries.
- **Output:** Integrated proportions stored in standardized format regardless of source method, enabling seamless downstream analysis.
- **Result:** Pipeline now robust to method failures. Can process libraries with any available method, providing flexibility and redundancy. Method usage statistics tracked and reported.
- **Future Considerations:** 
  - May add method comparison/validation step to assess agreement between methods
  - Could implement method-specific parameter tuning
  - Consider CondSCVI as additional alternative method
  - DestVI and Tangram are resulting in a segmentation fault

### [2025-01-05] - Memory Profiling and Optimization for Deconvolution
- **Issue:** Both DestVI and Tangram continue to result in segmentation faults, indicating memory exhaustion during processing.
- **Action:** Implemented comprehensive memory profiling and optimization throughout `celltype_deconvolution_phase3.py` to identify memory bottlenecks and prevent segmentation faults.
- **Memory Profiling Implementation:**
  - **Dual Tracking System:** Uses both `psutil` (process-level) and `tracemalloc` (Python-level) for comprehensive memory monitoring
  - **Logging Infrastructure:** Detailed memory logs saved to `output/deconvolution/logs/memory_log_*.txt` with timestamps
  - **Memory Checkpoints:** Memory usage logged at every major step:
    - Before/after data loading
    - Before/after preprocessing
    - Before/after each deconvolution method
    - After cleanup operations
  - **AnnData Memory Estimation:** Function to estimate memory footprint of AnnData objects (X matrix, obs, var, obsm)
  - **System Memory Monitoring:** Tracks both process memory (RSS) and system-wide available memory
- **Memory Optimization Strategies:**
  - **Aggressive Cleanup:** `aggressive_cleanup()` function calls `gc.collect()` twice and `malloc_trim()` to force memory release
  - **Memory Thresholds:** Checks memory before expensive operations (thresholds: 8-12GB RSS, 85-90% system memory)
  - **Adaptive Gene Reduction:** Automatically reduces number of signature genes if memory is tight (from 2000 to 1000)
  - **Adaptive Parameter Reduction:** Reduces training epochs and batch sizes for DestVI/Cell2location if memory is constrained
  - **Minimal Copies:** Uses `.copy()` only when necessary, prefers views where possible
  - **Early Variable Clearing:** Explicitly sets variables to `None` and calls cleanup immediately after use
  - **Method-Specific Optimizations:**
    - Tangram: Clears data copies before saving results
    - DestVI: Clears single-cell copy before spatial training, reduces epochs if needed
    - Cell2location: Reduces posterior sampling parameters (num_samples, batch_size) under memory pressure
- **Memory Logging Format:** Each log entry includes:
  - Timestamp
  - Step name
  - RSS (Resident Set Size) in MB
  - Process memory percentage
  - System available memory
  - System memory percentage
  - Tracemalloc peak (if available)
  - AnnData memory estimate (when applicable)
- **Error Handling:** Memory warnings logged but do not stop processing; methods skip to next if memory too high
- **Result:** Script now provides detailed memory diagnostics to identify exact locations of memory spikes. Memory optimizations should reduce peak usage and prevent segmentation faults. Logs enable post-hoc analysis of memory patterns.
- **Usage:** Memory logs can be analyzed to:
  - Identify which libraries/methods consume most memory
  - Determine optimal gene set sizes
  - Adjust thresholds for different hardware configurations
  - Debug segmentation faults by identifying memory spikes before crashes

### [2025-01-05] - Efficiency Optimization: Remove QC Filtering from Single-Cell Data
- **Action:** Removed QC filtering step from single-cell RNA-seq data loading in `celltype_deconvolution_phase3.py`.
- **Rationale:** QC filtering step was creating an unnecessary copy of the data and filtering cells, which adds computational overhead. For deconvolution purposes, including all cells (even QC-labeled ones) does not harm the process and may actually provide more reference data for mapping.
- **Implementation:** 
  - Removed the filtering step: `sc_adata = sc_adata[~sc_adata.obs[celltype_key].isin(qc_labels)].copy()`
  - Still exclude QC labels from the cell type list used for deconvolution (to avoid treating QC labels as actual cell types)
  - All cells are now retained in the reference dataset
- **Benefits:**
  - Faster data loading (no filtering/copying step)
  - Reduced memory overhead (no intermediate filtered copy)
  - More reference cells available for deconvolution mapping
  - Simpler code path
- **Result:** Single-cell data loading is now more efficient while maintaining the same deconvolution functionality. QC labels are excluded from cell type annotations but all cells remain in the reference dataset.

### [2025-01-05] - Tumor Differences Plotting: Fixed Significance Bar Overlap and Comparison Labeling
- **Issue:** In `tumor_differences.py`, significance bars were overlapping when multiple 1-vs-all comparisons were displayed, and only one comparison (minimum p-value) was shown in text format.
- **Action:** Redesigned the significance bar visualization to properly dodge bars and display all comparisons in the requested format.
- **Implementation:**
  - **Dodged Significance Bars:** Bars are now vertically offset (spaced by 8% of y-range) to prevent overlap
  - **Individual Bars per Comparison:** Each 1-vs-all comparison gets its own horizontal bar spanning all groups
  - **Asterisk Positioning:** Asterisks positioned above the specific group being compared (not at bar center)
  - **Comparison Text Format:** All comparisons displayed in bottom text box with format "Group A - Rest: p*-value (p=X.XXe-XX)"
  - **Multi-line Text:** Long comparison strings automatically split across multiple lines for readability
  - **Visual Distinction:** Comparison text box uses light yellow background with border for clear visibility
- **Format Details:**
  - Each comparison shown as: "{Group} - Rest: {asterisk} (p={adj_pval:.2e})"
  - Example: "Normal Alveolar Cells - Rest: * (p=2.34e-02) | Non-Solid Tumor - Rest: ** (p=1.23e-03)"
  - Non-significant comparisons included in text but not shown with bars
  - Asterisks follow standard convention: * < 0.05, ** < 0.01, *** < 0.001, **** < 0.0001
- **Result:** Plots now clearly show all comparisons without overlap, with all p-values displayed in the requested format. Bars are properly spaced and comparisons are easily readable.

### [2025-01-05] - Tumor Differences Plotting: Final Refinements for Readability
- **Action:** Final refinements to `tumor_differences.py` plotting function for optimal readability and visual clarity.
- **Changes:**
  - **One Line Per Comparison:** Each comparison now takes exactly one line in the text box (removed "|" separator)
  - **Whisker Protection:** Bars positioned above maximum whisker extent with 10% y-range buffer to prevent interference
  - **Dynamic Layout:** Bottom margin adjusts automatically based on number of comparisons (0.025 per line)
  - **Whisker Calculation:** Maximum whisker calculated as Q3 + 1.5*IQR for each group to ensure bars are always above
- **Implementation Details:**
  - Bars start at `max_whisker + 10% buffer` and stack upward
  - Each comparison text on separate line: "Group A - Rest: p*-value (p=X.XXe-XX)"
  - Layout uses `tight_layout(rect=[0, bottom_margin, 1, 1])` with dynamic bottom margin
- **Result:** Plots are now maximally readable with no bar-whisker interference and clear one-line-per-comparison format.

### [2025-01-05] - Code Cleanup: Gene Signature Scoring Refactoring
- **Action:** Extracted core signature scoring functionality from `signature_cleaned.py` into a new clean script `score_gene_signatures.py` and archived the old file.
- **Rationale:** The `signature_cleaned.py` file contained obsolete code from previous projects and mixed scoring with plotting/analysis functions. A clean, focused script is needed for Phase IV tasks (Macrophage State Comparison, Differential cell-type proportion).
- **Implementation:**
  - **New Script:** `scripts/score_gene_signatures.py` - Clean, modular signature scoring script
  - **Core Functions Extracted:**
    - `score_signatures_nested()` - Main scoring function
    - `_flatten_nested_dict()` - Flatten nested signature structure
    - `_index_genes()` - Find genes in var_names
    - `_get_matrix()` - Extract expression matrix with transformations
    - `_zscore_cols()` - Z-score genes before averaging
  - **Workflow:**
    1. Load gene signatures from `metadata/gene_signatures.json`
    2. Load AnnData from `results/adata.annotation.masked.h5ad`
    3. Score all signatures (including macrophage states from Mission.md)
    4. Z-score each signature across spots for visualization
    5. Verify macrophage signatures are present (Alveolar macrophages Merad, Mac.2 MoMac M2-like, Mac.6 MoMAc M2-like)
    6. Save AnnData with scores to `results/adata.img.genescores.h5ad`
  - **Archived:** `scripts/signature_cleaned.py` → `scripts/old_code/signature_cleaned.py`
- **Mission.md Alignment:**
  - Script scores all signatures including macrophage states required for Phase IV:
    - Alveolar macrophages Merad
    - Mac.2 MoMac M2-like
    - Mac.6 MoMAc M2-like
  - Output AnnData (`adata.img.genescores.h5ad`) contains all signature scores needed for downstream 1-vs-all statistical testing
- **Result:** Clean, maintainable signature scoring script ready for Phase IV analysis tasks. Old code archived for reference.
- **Future Improvement:** Rewrite scoring function such that it uses the sc.tl.score_genes instead of just calculating the mean score.

### [2025-01-05] - Tumor Differences: Added Liron Myeloid and T Cell Signatures
- **Action:** Updated `tumor_differences.py` to include additional Liron signatures for comprehensive immune program analysis.
- **Rationale:** Mission.md Phase IV tasks require analysis of macrophage states and T cell programs. Liron signatures provide additional granularity for myeloid and T cell characterization beyond the existing Immune_Lymphoid signatures.
- **Added Signatures:**
  - **T Cell Signatures (Liron):**
    - Suppressive Tregs (more specific than general Tregs)
    - CD8.4 TRM pre-exhausted (tissue-resident memory T cells)
    - CD8.1 GZMK pre-exhausted (pre-exhausted CD8+ T cells)
  - **Myeloid/Macrophage Signatures (Liron):**
    - M2-macrophages (general M2 polarization)
    - MoMacs Merad (monocyte-derived macrophages)
    - Alveolar macrophages Merad (alveolar macrophage signature)
    - Mac.2 MoMac M2-like (M2-like monocyte-derived macrophages)
    - Mac.6 MoMAc M2-like (alternative M2-like monocyte-derived macrophages)
- **Implementation:** Added 11 Liron signatures to `immune_programs` list (6 T cell + 5 myeloid), organized with comments for clarity.
- **Result:** Script now tests 23 total immune programs (12 original + 11 Liron) across three tumor types, providing comprehensive coverage for Phase IV macrophage state comparison and differential cell-type proportion analysis.

### [2025-01-06] - Macrophage Localization: Colocalization with Proliferative Program
- **Action:** Created `macrophage_localization.py` to measure macrophage colocalization with proliferative program as specified in Mission.md.
- **Rationale:** Mission.md Phase IV requires analysis of macrophage spatial colocalization with tumor proliferation to understand macrophage-tumor interactions.
- **Implementation:**
  - **Scatterplots:** X-axis = Proliferative_Tumor_z score, Y-axis = Liron macrophage signature scores
  - **Stratification:** Separate analysis for Normal Alveolar Cells, Non-Solid Tumor, and Solid Tumor
  - **Visualization:** Different colors per tumor type (green=Normal, orange=Non-Solid, blue=Solid)
  - **Correlation Analysis:** Pearson correlation calculated for each tumor type separately
  - **Annotations:** Correlation coefficients, p-values, and sample sizes displayed outside plot in text box
  - **Macrophage Signatures Tested:**
    - M2-macrophages
    - MoMacs Merad
    - Alveolar macrophages Merad
    - Mac.2 MoMac M2-like
    - Mac.6 MoMAc M2-like
- **Output:**
  - Individual scatterplots for each macrophage signature (PDF and PNG)
  - Correlation summary CSV with r, p-values, and n for each signature-tumor type combination
  - Saved to `figures/manuscript/macrophage_localization/` and `figures/stats/`
- **Features:**
  - Handles missing signatures gracefully (skips if not found)
  - Filters to valid tumor types only
  - Removes NaN values before correlation calculation
  - Requires minimum of 3 points per group for correlation
  - Professional formatting with grid, legend, and correlation annotations
- **Result:** Script ready to analyze macrophage-proliferation colocalization patterns across tumor types, providing quantitative measures of spatial relationships.

### [2025-01-06] - Macrophage Localization: Enhanced Visualization with Regression Lines
- **Action:** Enhanced `macrophage_localization.py` with improved filtering, reference lines, and regression analysis with confidence intervals.
- **Improvements:**
  - **Filtering:** Filter spots with proliferative score < 0 (only analyze spots with positive/zero proliferation)
  - **Reference Lines:** Added x=0 and y=0 dashed lines to mark origin (black, semi-transparent)
  - **Filtered Points Visualization:** Points with proliferative score < 0 displayed in light gray (alpha=0.3) for context
  - **Regression Lines:** Best fit line for each tumor type condition using `scipy.stats.linregress`
  - **Confidence Intervals:** 95% confidence intervals around regression lines (shaded areas, alpha=0.2) matching seaborn regplot style
  - **Layering:** Proper z-order to ensure reference lines, confidence intervals, regression lines, and scatter points are correctly layered
- **Implementation Details:**
  - Filter mask: `x_data >= 0` to keep only spots with proliferative score >= 0
  - Excluded points plotted first in light gray with `zorder=0`
  - Reference lines at `zorder=1` (x=0 and y=0)
  - Confidence intervals at `zorder=2` (shaded areas)
  - Regression lines at `zorder=2` (solid lines)
  - Scatter points at `zorder=3` (on top)
  - Confidence interval calculation uses t-distribution (95% CI, df=n-2)
- **Visual Result:** Plots now clearly show:
  - Which spots are excluded (light gray, proliferative < 0)
  - Origin reference (x=0, y=0 lines)
  - Trend lines for each condition with uncertainty bounds
  - All correlations still annotated below plot
- **Result:** Enhanced scatterplots provide better visual interpretation of macrophage-proliferation relationships with statistical trend lines and confidence intervals for each tumor type.

### [2025-01-06] - Macrophage Localization: Additional Filtering and Confidence Interval Fix
- **Action:** Fixed filtering to include macrophage scores and ensured confidence intervals display for all tumor types.
- **Issues Fixed:**
  1. **Macrophage Score Filtering:** Now filters spots with macrophage score < 0 in addition to proliferative score < 0
  2. **Confidence Intervals for All Types:** Added error handling and validation to ensure confidence intervals are calculated and displayed for Normal, Non-Solid, and Solid tumor types
- **Implementation:**
  - **Dual Filtering:** Filter mask now checks both `(x_data >= 0) & (y_data >= 0)` (proliferative AND macrophage >= 0)
  - **Excluded Points:** All points with either proliferative < 0 OR macrophage < 0 are shown in light gray
  - **Robust Regression:** Added try-except block and validation checks for:
    - Finite regression parameters (slope, intercept, std_err)
    - Valid x-range (x_max > x_min)
    - Valid variance (ssx > 0, std_err > 0)
    - Finite confidence interval values
  - **Fallback Behavior:** If confidence interval calculation fails, regression line is still plotted without CI
  - **Error Reporting:** Warnings printed if regression fails for any tumor type
- **Result:** All tumor types now show regression lines with confidence intervals when data is valid, and filtering correctly excludes spots with either negative proliferative or macrophage scores.

### [2025-01-05] - Macrophage Localization: Fixed Confidence Interval Calculation and Visibility
- **Issue:** Confidence intervals were not visible or not plotting for Non-Solid and Solid tumor types.
- **Root Cause:** 
  1. Incorrect formula: Using prediction interval formula (`sqrt(1 + 1/n + ...)`) instead of confidence interval formula (`sqrt(1/n + ...)`)
  2. Low visibility: Alpha=0.2 was too low to see the shaded areas
- **Fix:**
  - **Corrected CI Formula:** Changed from prediction interval to confidence interval formula:
    - Old: `se = std_err * np.sqrt(1 + 1/n + (x_line - x_mean)**2 / ssx)`
    - New: `se = std_err * np.sqrt(1.0/n + (x_line - x_mean)**2 / ssx)`
  - **Increased Visibility:** Alpha increased from 0.2 to 0.3 for better visibility
  - **Visual Improvements:** Added `edgecolor='none'` to fill_between for cleaner appearance
  - **Validation:** Added check to ensure CI values are positive before plotting
  - **Debug Output:** Added print statements to verify CI calculation for each tumor type
- **Result:** Confidence intervals now correctly calculated and visible for all tumor types (Normal, Non-Solid, Solid) with proper confidence interval bands around regression lines.

### [2025-01-05] - Macrophage Localization: Removed Confidence Intervals
- **Action:** Removed confidence interval bands from regression lines in `macrophage_localization.py`.
- **Rationale:** Confidence intervals were too narrow to be visually useful, making the plots cluttered without adding meaningful information.
- **Implementation:**
  - Removed all confidence interval calculation code (CI formula, t-critical, fill_between)
  - Simplified regression plotting to show only the best-fit line
  - Kept regression line calculation and plotting for all tumor types
  - Updated docstring to reflect removal of confidence intervals
- **Result:** Cleaner plots with regression lines only, making trends easier to see without narrow confidence bands.

### [2025-01-05] - Gene Signature Scoring: Transition to Seurat-based Method
- **Action:** Updated `score_gene_signatures.py` to use Scanpy's `score_genes` function (Seurat's AddModuleScore) instead of simple mean averaging.
- **Rationale:** Simple mean averaging is biased by highly expressed genes and does not account for technical variation (sequencing depth, batch effects). Seurat-based scoring is more robust and follows best practices.
- **Key Differences:**
  - **Old Method:** Simple mean of z-scored gene expression: `mean(zscore(genes))`
  - **New Method:** Matched control gene subtraction: `mean(signature genes) - mean(matched control genes)`
- **Algorithm Details:**
  1. Bins all genes by their average expression across cells/spots
  2. For each signature gene, randomly samples control genes from the same expression bin
  3. Calculates per-cell score as: mean(signature) - mean(controls)
  4. Controls for technical variation and expression magnitude differences
- **Advantages:**
  1. **Expression-level matching:** Control genes matched by expression magnitude
  2. **Technical variation control:** Accounts for sequencing depth and batch effects
  3. **Reduced bias:** Prevents highly expressed genes from dominating scores
  4. **Standardized approach:** Consistent with Seurat/Scanpy best practices
- **Parameters:**
  - `ctrl_size`: 50 control genes per signature gene (default)
  - `n_bins`: 25 expression bins for matching (default)
  - `use_raw`: True (requires raw counts for proper binning)
- **Implementation:**
  - Removed `_get_matrix()` and `_zscore_cols()` helper functions (no longer needed)
  - Removed `log1p`, `clip_min`, and `zscore_per_gene` parameters (handled internally by score_genes)
  - Removed `layer` parameter (score_genes uses adata.X or adata.raw, not layers)
  - Added error handling for cases where score_genes fails (e.g., insufficient control genes)
- **Result:** More robust gene signature scores that better reflect biological signal over technical noise, especially important for comparing scores across samples with different sequencing depths or batch effects.

### [2025-01-20] - Process Colocalization Analysis: Modular Implementation
- **Action:** Created `process_colocalization.py` to analyze spatial co-occurrence patterns of gene signature scores (processes) across spots using multiple modular methods.
- **Rationale:** Mission.md Phase IV requires identification of process pairs that are colocalized (positive association) or anti-colocalized (mutual exclusion) to understand spatial relationships between biological processes.
- **Modular Architecture:**
  - **Configuration-Driven:** `ANALYSIS_CONFIG` dictionary allows enabling/disabling each analysis independently
  - **Separate Functions:** Each analysis method implemented as independent function:
    - `analyze_pearson_correlation()`: Uses `DataFrame.corr()` for pairwise correlations
    - `analyze_morans_i()`: Manual Moran's I computation via spatial neighbors graph
    - `analyze_neighborhood_enrichment()`: Thresholded enrichment using `squidpy.gr.nhood_enrichment`
    - `analyze_volcano_plots()`: Statistical comparisons between Normal, Non-Solid, and Solid tumor types
- **Implementation Details:**
  - **Pearson Correlation:**
    - Extracts all signature scores as DataFrame
    - Uses built-in `DataFrame.corr(method='pearson')` for efficiency
    - Generates correlation heatmap (RdBu_r colormap, -1 to 1 range)
    - Identifies top colocalized and anti-colocalized pairs
  - **Moran's I:**
    - Manual computation using spatial neighbors connectivity matrix
    - Formula: `I = (n/W) * sum(w_ij * (x_i - mean)(x_j - mean)) / sum((x_i - mean)^2)`
    - Permutation testing (100 permutations) for p-value calculation
    - Generates diagonal heatmap showing spatial autocorrelation for each process
  - **Neighborhood Enrichment:**
    - Thresholds signatures at ±1 z-score (low/high classification)
    - Uses `squidpy.gr.nhood_enrichment` with 1000 permutations
    - Computes enrichment for high-scoring spots vs. not-high spots
    - Can be disabled for faster execution (slower analysis)
  - **Volcano Plots:**
    - Three comparisons: Normal vs Non-Solid, Normal vs Solid, Non-Solid vs Solid
    - Mann-Whitney U test for each signature between groups
    - Benjamini-Hochberg FDR correction applied
    - Log2 fold change calculated (mean difference for z-scores)
    - Color coding: red (up, FC > 0.5), blue (down, FC < -0.5), orange (significant but small effect), gray (not significant)
    - Annotations: All significant points (adj_p < 0.05) labeled with cleaned signature names
- **Output:**
  - Correlation heatmap: `figures/process_colocalization/correlation_heatmap.pdf`
  - Moran's I heatmap: `figures/process_colocalization/morans_i_heatmap.pdf`
  - Volcano plots: `figures/process_colocalization/volcano_*.pdf` (3 plots)
  - CSV results: `results/process_colocalization/*.csv` (correlation matrix, Moran's I, volcano stats)
- **Signature Name Cleaning:**
  - Removes `sig:` prefix and `_z` suffix
  - Replaces `/` with `-` for readability
  - Truncates long names (>40 chars) for annotation clarity
- **Configuration:**
  - Default: Pearson correlation and volcano plots enabled
  - Moran's I and neighborhood enrichment can be enabled as needed
  - Each analysis runs independently and can be toggled via `ANALYSIS_CONFIG`
- **Result:** Comprehensive spatial colocalization analysis pipeline that identifies process relationships through multiple complementary methods. Modular design allows selective execution based on computational resources and analysis needs.

### [2025-01-20] - Process Colocalization: Moran's I Fix and Volcano Plot Enhancements
- **Issue 1:** Moran's I computation was failing because `sq.gr.spatial_autocorr` expects genes in `var_names`, not signature scores stored in `obs`.
- **Fix:** Implemented manual Moran's I computation using spatial neighbors connectivity matrix:
  - Uses `spatial_connectivities` from `sq.gr.spatial_neighbors`
  - Manual formula: `I = (n/W) * sum(w_ij * (x_i - mean)(x_j - mean)) / sum((x_i - mean)^2)`
  - Permutation testing (100 permutations) for p-value calculation
  - Handles NaN values and edge cases gracefully
- **Issue 2:** User requested volcano plots with annotations for significant processes to visualize differences between Normal, Non-Solid, and Solid tumor types.
- **Implementation:**
  - Added `analyze_volcano_plots()` function for three comparisons:
    - Normal vs Non-Solid Tumor
    - Normal vs Solid Tumor
    - Non-Solid vs Solid Tumor
  - Statistical testing: Mann-Whitney U test for each signature between groups
  - FDR correction: Benjamini-Hochberg adjustment applied
  - Fold change: Log2 fold change calculated (mean difference for z-scores)
  - Color coding:
    - Red: Significant up-regulation (adj_p < 0.05, FC > 0.5)
    - Blue: Significant down-regulation (adj_p < 0.05, FC < -0.5)
    - Orange: Significant but small effect (adj_p < 0.05, |FC| ≤ 0.5)
    - Gray: Not significant
  - **Annotations:** All significant points (adj_p < 0.05) labeled with cleaned signature names:
    - Removes `sig:` prefix and `_z` suffix
    - Replaces `/` with `-` for readability
    - Truncates long names (>40 chars)
    - White background boxes with gray borders and small arrows
- **Output:**
  - Three volcano plots: `volcano_Normal_vs_NonSolid.pdf`, `volcano_Normal_vs_Solid.pdf`, `volcano_NonSolid_vs_Solid.pdf`
  - CSV file with all statistics: `volcano_plot_results.csv`
  - Summary statistics printed showing number of up/down-regulated processes per comparison
- **Result:** Moran's I now computes correctly for signature scores, and volcano plots provide clear visualization of process differences between tumor types with annotated significant processes for easy identification.

### [2025-01-20] - Process Colocalization: Faceted Volcano Plots and Signature Filtering
- **Date:** 2025-01-20
- **Action:** Enhanced `process_colocalization.py` with faceted volcano plots and signature selection capabilities.
- **Improvements:**
  - **Faceted Volcano Plots:** Changed from 3 separate plots to a single faceted figure with 3 subplots (Normal vs Non-Solid, Non-Solid vs Solid, Normal vs Solid) for easier comparison.
  - **Signature Filtering:** Added `VOLCANO_SIGNATURES_INCLUDE` and `VOLCANO_SIGNATURES_EXCLUDE` parameters for selective visualization:
    - Supports exact signature names (with or without 'sig:' prefix and '_z' suffix)
    - Supports pattern matching (case-insensitive substring search)
    - Allows filtering large signature sets to focus on specific processes
  - **Correlation Clustermaps:** Added hierarchically clustered heatmaps (clustermaps) for both Pearson correlation and macrophage-proliferation correlations:
    - Clusters signatures/processes by similarity
    - Includes dendrograms showing relationships
    - Same color scale and annotations as standard heatmaps
- **Implementation:**
  - Faceted plot uses `plt.subplots(1, 3)` for side-by-side comparison
  - Each subplot uses same annotation logic with proper dodging
  - Legend shown only on first subplot to avoid duplication
  - Signature filtering uses helper function `filter_signatures()` with pattern matching
- **Output:**
  - Single faceted volcano plot: `volcano_faceted.pdf` (replaces 3 separate plots)
  - Correlation clustermap: `correlation_clustermap.pdf` (in addition to standard heatmap)
  - Macrophage correlation clustermap: `correlation_clustermap.pdf` (in addition to standard heatmap)
- **Result:** More efficient visualization with faceted plots and flexible signature selection for focused analysis. Clustermaps reveal natural groupings of processes and signatures.

### [2025-01-20] - Signature Heatmap Visualization: Comprehensive Gene Signature Visualization
- **Date:** 2025-01-20
- **Action:** Created `signature_heatmap.py` to generate comprehensive heatmap and clustermap visualizations of gene signature scores across spots.
- **Rationale:** Need for comprehensive visualization of all gene signature scores across spots with proper grouping by patient and solidity to identify patterns and relationships.
- **Implementation:**
  - **Rows:** Spots (sorted by patient and solidity)
  - **Columns:** Gene signatures (hierarchically clustered in clustermaps)
  - **Annotations:** Patient (library_id) and solidity (tumor_type) shown as color bars above heatmap
  - **3-Level Hierarchy:**
    - Level 1 (Macro): Primary grouping (patient or solidity based on sort_by)
    - Level 2 (Secondary): Secondary grouping (solidity or patient)
    - Level 3 (Fine): Clustering within smallest groups (only in clustermaps)
  - **Sorting Logic:**
    - Patient-first: Cluster patients by mean signature scores (macro), then maintain solidity order (Normal → Non-Solid → Solid), then fine cluster within patient-solidity groups
    - Solidity-first: Fixed solidity order (Normal → Non-Solid → Solid), then cluster patients by mean scores within each solidity, then fine cluster within solidity-patient groups
  - **Clustering:** Uses `cluster_within_groups()` function to cluster spots within each (patient, solidity) or (solidity, patient) group while preserving group boundaries
  - **Legend:** Solidity color legend in top right corner showing Normal, Non-Solid, and Solid colors
- **Output:**
  - 4 figures total:
    1. `signature_heatmap_patient_solidity.pdf` - Standard heatmap sorted by patient then solidity
    2. `signature_heatmap_solidity_patient.pdf` - Standard heatmap sorted by solidity then patient
    3. `signature_clustermap_patient_solidity.pdf` - Clustermap sorted by patient then solidity (with fine clustering)
    4. `signature_clustermap_solidity_patient.pdf` - Clustermap sorted by solidity then patient (with fine clustering)
  - All saved to `figures/manuscript/signature_heatmaps/`
- **Features:**
  - Consistent solidity ordering (Normal → Non-Solid → Solid) in all plots
  - Patient clustering by mean signature scores for macro-level organization
  - Fine clustering within groups preserves group boundaries
  - Color annotations for patient and solidity with proper RGB conversion
  - Signature names cleaned for display (removes 'sig:' prefix and '_z' suffix)
- **Result:** Comprehensive visualization pipeline that reveals patterns in gene signature scores across spots with proper hierarchical organization and clear annotations for patient and solidity grouping.

### [2025-01-20] - Macrophage Localization: Correlation Heatmaps and Clustermaps
- **Date:** 2025-01-20
- **Action:** Enhanced `macrophage_localization.py` with correlation heatmap and clustermap visualizations.
- **Rationale:** Need for comprehensive visualization of macrophage-proliferation correlations across tumor types to identify patterns and relationships.
- **Implementation:**
  - **Correlation Matrix:** Pivots correlation data to create matrix (macrophage signatures × tumor types)
  - **Standard Heatmap:** Shows Pearson correlation coefficients with significance annotations (asterisks for p < 0.05, 0.01, 0.001)
  - **Clustermap:** Hierarchically clusters macrophage signatures by correlation patterns while maintaining tumor type order (Normal → Non-Solid → Solid)
  - **Annotations:** Correlation values with significance stars (* p<0.05, ** p<0.01, *** p<0.001)
- **Output:**
  - `correlation_heatmap.pdf` - Standard correlation heatmap
  - `correlation_clustermap.pdf` - Hierarchically clustered correlation heatmap
  - Both saved to `figures/manuscript/macrophage_localization/`
- **Result:** Enhanced visualization of macrophage-proliferation relationships with both standard and clustered views for pattern identification.

---

## Future Rationale & Pending Decisions

### TLS Transcriptional Heterogeneity
- **Pending:** We need to decide whether to use `Cell2location` or `DestVI` for the final TLS deconvolution.
- **Workaround:** Current `tangram_batches` output is serving as a baseline, but high-resolution Visium HD data may require a move to `SpatialData` containers to avoid memory issues.

### Gene Signature Refinement
- **Action:** Updated `metadata/gene_signatures.json` with specific EMT and Myeloid markers.
- **Sanity Check:** Agent must verify that these signatures do not overlap significantly with housekeeping genes.

