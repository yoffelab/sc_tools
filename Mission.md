# Mission: Lung Tumor Evolution and TLS Transcriptomics

**Current Status:** Phase III-IV (Cell Type Deconvolution, Spatial Structures & Niche Analysis)
**Author:** Junbum Kim
**Last Updated:** 2025-01-20

---

## 1. Objective
Identify transcriptional and spatial differences across the lung tumor progression spectrum (Normal, Non-Solid, and Solid tumors) with a specific focus on the heterogeneity of Tertiary Lymphoid Structures (TLS).

---

## 2. Completed Tasks
- [x] **Data Ingestion:** Converted raw Visium/scRNA-seq data into AnnData containers (`results/adata_filtered.h5ad`).
- [x] **Latent Modeling:** Trained scVI models to generate integrated latent spaces for batch-corrected analysis (`results/scvi.h5ad`).
- [x] **Gene Scoring:** Calculated LUAD-related scores using Seurat-based method (Scanpy `score_genes`) from established signatures (`results/adata.img.genescores.h5ad`, using `metadata/gene_signatures.json`).
- [x] **Initial Visualization:** Generated per-sample QC and EMT spatial distribution plots.
- [x] **Grouping:** Categorized samples by pathology (Normal, Non-Solid, Solid) in `adata.obs`.
- [x] **Differential Program Analysis:** Performed 1-vs-all statistical testing on tumor and immune programs across tumor types with boxplots and significance annotations (`scripts/tumor_differences.py`).
- [x] **Macrophage Localization:** Analyzed macrophage colocalization with proliferative program, stratified by tumor type (`scripts/macrophage_localization.py`).
- [x] **Spot Process Colocalization Analysis:** Implemented modular analysis pipeline for spatial co-occurrence patterns of processes (`scripts/process_colocalization.py`).
- [x] **Signature Heatmaps:** Created comprehensive heatmap and clustermap visualizations of gene signature scores with patient and solidity annotations (`scripts/signature_heatmap.py`).
---

## 3. Active Tasks (Roadmap)
### Phase III: Cell Type Deconvolution
- [x] **Deconvolution Pipeline:** Implemented modular deconvolution pipeline with memory optimization (`scripts/celltype_deconvolution_phase3.py`). *critical* - still need to identify minimum required unit for the method to work, and the size of data that causes things to fail. 
    - [x] `tangram`: Implemented with batch processing per library_id to handle memory constraints.
    - [ ] `cell2location`: Test implementation with batch processing to identify memory breaking points.
    - [ ] `DestVI`: Test implementation with batch processing.
    - [x] **Fallback Strategy:** Implemented method fallback (DestVI → Cell2location → Tangram) with error handling.

- [ ] **Gene Signature Refinement:** Update `metadata/gene_signatures.json` to contain specific and accurate gene lists for associated processes.
    - [ ] Add relevant processes (e.g., Hallmark gene sets).
    - [ ] Update gene lists to be robust, specific, and accurate.
    - [ ] Validate against external databases (HLCA, MSigDB, TCGA LUAD).

### Phase IV: Comparative Tumor Analysis
- [x] **Differential Program Analysis:** (`scripts/tumor_differences.py`)
    - Performed 1-vs-all statistical testing on tumor programs (EMT, Hypoxia, Proliferative) across tumor types.
    - Performed 1-vs-all statistical testing on immune programs (Immune_Lymphoid, Innate_Other, Processes:TLS_Formation, Liron signatures) across tumor types.
    - Generated boxplots with significance asterisks, adjusted p-values, and comparison text annotations.
- [ ] **Spatially Variable Genes (SVG):** (`scripts/spatial_analysis.py`)
    - Identify genes that characterize the transition from Non-Solid to Solid tumors using Moran's I or `squidpy.gr.spatial_autocorr`.
    - Generate spatial plots highlighting SVG patterns.
- [x] **Spot Process Colocalization Analysis:** (`scripts/process_colocalization.py`)
    - **Objective:** Analyze spatial co-occurrence patterns of processes across spots.
    - **Methods Implemented:**
        - [x] Pearson correlation of processes across spots (using `DataFrame.corr()`).
        - [x] Moran's I for spatial autocorrelation of processes (manual computation via spatial neighbors graph).
        - [ ] Thresholded neighborhood enrichment using `squidpy.gr.nhood_enrichment` after thresholding processes by low / high (-1, +1) (implemented, disabled by default).
        - [x] Volcano plots for Normal vs Non-Solid vs Solid comparisons with statistical testing and annotations (faceted plot with 3 subplots).
        - [x] Signature filtering: Include/exclude parameters for selecting which signatures to plot.
    - **Output:** 
        - Correlation heatmap and clustermap identifying colocalized (positive) and anti-colocalized (negative) process pairs.
        - Moran's I heatmap showing spatial autocorrelation for each process.
        - Faceted volcano plots with annotated significant processes for each tumor type comparison.
        - CSV files with all statistical results.
- [x] **Signature Heatmap Visualization:** (`scripts/signature_heatmap.py`)
    - **Objective:** Create comprehensive heatmap and clustermap visualizations of gene signature scores across spots.
    - **Features:**
        - Rows: Spots (sorted by patient and solidity)
        - Columns: Gene signatures (hierarchically clustered in clustermaps)
        - Annotations: Patient (library_id) and solidity (tumor_type) color bars
        - 3-level hierarchy: Macro grouping (patient/solidity), secondary grouping, fine clustering within groups
    - **Output:**
        - 4 figures total: 2 heatmaps (patient→solidity, solidity→patient) and 2 clustermaps (same sorting)
        - Solidity legend in top right corner
        - All figures saved to `figures/manuscript/signature_heatmaps/`
- [ ] **Cell Type Colocalization Analysis:**
    - **Objective:** Analyze spatial co-occurrence patterns of cell type proportions (from deconvolution) across spots.
    - **Methods:**
        - [ ] Pearson correlation of cell type proportions across spots.
        - [ ] Moran's I for spatial autocorrelation of cell type proportions.
        - [ ] Thresholded neighborhood enrichment using `squidpy.gr.nhood_enrichment`.
        - [ ] Tensor decomposition (MEFISTO) for spatiotemporal patterns (if applicable).
    - **Output:** Identify cell type pairs that are colocalized (positive association) or anti-colocalized (mutual exclusion).
- [ ] **Spatial Gene Expression Patterns:**
    - **Objective:** Identify gene signatures that up/down-regulate as a result of spatial proximity between cell types.
    - **Approach:** Explore graph-based methods where:
        - Nodes represent spots with cell type proportions.
        - Edges represent spatial proximity with associated gene signature scores.
    - **Implementation:** Research existing methods (e.g., spatial graph neural networks) or develop algorithm for expected/observed transcriptional programs.
- [ ] **Spatial Transition Areas:**
    - Identify transcriptional upregulation/downregulation in transition areas:
        - Between Normal, Non-Solid, and Solid tumor regions.
        - Between tumor cores and TLS.
    - Use spatial distance metrics and differential expression analysis.
- [ ] **Macrophage State Comparison:**
    - Compare macrophage subpopulations across tumor types:
        - Alveolar macrophages Merad
        - Monocyte-derived macrophages (M1/M2 polarization)
        - M2 subpopulations: `Mac.2 MoMac M2-like`, `Mac.6 MoMAc M2-like`
    - Perform 1-vs-all statistical testing on macrophage signatures across tumor types.
    - Generate boxplots with significance annotations.
- [ ] **Differential Cell Type Proportions:**
    - Perform 1-vs-all statistical testing on cell type proportions (from deconvolution) across tumor types.
    - Generate boxplots with significance asterisks and adjusted p-values.
    
### Phase V: TLS-Specific Transcriptomics
- [ ] **TLS Niche Extraction:** (`scripts/tls_analysis.py`)
    - Subset the `tls_clustered.h5ad` object to focus on lymphoid-rich neighborhoods.
    - Generate spatial plots of TLS distribution across samples.
- [x] **B-cell/T-cell State Comparison:** (`scripts/tls_analysis.py`)
    - Analyzed ratios of Naive vs. Memory B-cells and Exhausted vs. Effector T-cells within TLS across tumor stages.
    - Performed 1-vs-all statistical testing with significance annotations.
- [x] **Ligand-Receptor Crosstalk:** (`scripts/tls_analysis.py`)
    - Used `squidpy.gr.ligrec` to infer signaling between tumor-adjacent TLS and tumor cells in Non-Solid vs. Solid samples.
- [x] **Macrophage Localization:** (`scripts/macrophage_localization.py`)
    - Measured extent of macrophage colocalization with proliferative program.
    - Stratified by normal, non-solid, and solid tumors.
    - Generated scatterplots with tumor proliferation score (x-axis) vs. macrophage scores (y-axis).
    - Annotated Pearson correlation per condition.
    - Generated correlation heatmap and clustermap for macrophage-proliferation correlations across tumor types.
---

## 4. Blockers and Sanity Checks
- **Sanity Check:** Ensure that low-count Non-Solid (GGO) spots are not being over-filtered during comparative QC.
- **Statistical Correction:** All differential expression and multiple comparison results must use Benjamini-Hochberg (FDR) adjustment.
- **Memory Management:** Deconvolution scripts must use batch processing per library_id to prevent segmentation faults.

---

## 5. Technical Decision Log (Reference @Journal.md)
- **scVI Latent Dimensions:** Selected n=30 latent dimensions for integration (see Journal.md).
- **Comparison Mode:** Default to 1-vs-all (1-vs-rest) comparisons for tumor types to highlight stage-specific programs.
- **Gene Signature Scoring:** Transitioned from simple mean to Seurat-based method (Scanpy `score_genes`) for robust scoring that controls for technical variation.
- **Deconvolution Strategy:** Implemented modular pipeline with fallback (DestVI → Cell2location → Tangram) and batch processing per library_id for memory efficiency.
- **Process Colocalization:** Implemented modular analysis pipeline with Pearson correlation (DataFrame.corr), Moran's I (manual computation), and faceted volcano plots for identifying colocalized and anti-colocalized process pairs across tumor types. Includes signature filtering (include/exclude) for selective visualization.
- **Signature Heatmaps:** Created comprehensive visualization pipeline with 3-level hierarchical sorting (patient/solidity grouping, mean score clustering, fine clustering) for both heatmaps and clustermaps. Maintains consistent solidity ordering (Normal → Non-Solid → Solid) and patient clustering by mean signature scores.