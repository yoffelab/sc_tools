"""Phase 3 (preprocess): Normalize, benchmark integrations, and cluster robin cell data.

Reads results/adata.annotated.h5ad, runs the full integration benchmark workflow:
1. Subsample ~50k cells (stratified by sample)
2. Run all candidate methods (Harmony, BBKNN, ComBat, Scanorama, scVI, PCA baseline)
3. Score by batch correction quality
4. Apply the best method to the full 1.5M-cell dataset
5. Cluster and UMAP on the winning embedding
6. Generate post-integration QC report

Saves results/adata.normalized.h5ad and results/integration_method.txt.

Run from project root or via SLURM (GPU recommended for scVI).
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
REPO_ROOT = PROJECT_ROOT.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def main():
    import scanpy as sc

    from sc_tools.bm.integration import run_full_integration_workflow
    from sc_tools.pp.normalize import backup_raw, filter_genes_by_pattern, log_transform, normalize_total
    from sc_tools.pp.reduce import cluster, neighbors, pca, umap

    input_path = PROJECT_ROOT / "results" / "adata.annotated.h5ad"
    output_path = PROJECT_ROOT / "results" / "adata.normalized.h5ad"
    results_dir = PROJECT_ROOT / "results"

    if not input_path.exists():
        logger.error("adata.annotated.h5ad not found at %s", input_path)
        sys.exit(1)

    logger.info("Reading %s", input_path)
    adata = sc.read_h5ad(input_path)
    adata.obs_names_make_unique()
    logger.info("Loaded: %d cells x %d genes", adata.n_obs, adata.n_vars)

    # Filter MT/RP/HB genes before anything
    filter_genes_by_pattern(adata)
    logger.info("After gene pattern filter: %d genes", adata.n_vars)

    # HVG on raw counts (seurat_v3 requires raw counts)
    sc.pp.highly_variable_genes(
        adata, n_top_genes=3000, flavor="seurat_v3", batch_key="sample"
    )
    n_hvg = adata.var["highly_variable"].sum()
    logger.info("Selected %d HVGs (seurat_v3 on raw counts)", n_hvg)

    # Backup raw counts, then normalize
    backup_raw(adata)
    normalize_total(adata)
    log_transform(adata)
    logger.info("Normalized and log-transformed")

    # Subset to HVGs for integration
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    logger.info("HVG subset: %d cells x %d genes", adata_hvg.n_obs, adata_hvg.n_vars)

    # PCA on HVG subset
    pca(adata_hvg, n_comps=50)

    # Run integration benchmark: subsample 50k, test all methods, apply best
    logger.info("Starting integration benchmark workflow (subsample=50000)...")
    adata_hvg, comparison_df, best_method = run_full_integration_workflow(
        adata_hvg,
        modality="xenium",  # transcriptomic methods
        batch_key="sample",
        celltype_key=None,  # no celltypes yet
        methods=["harmony", "bbknn", "combat", "scanorama", "scvi", "pca"],
        output_dir=str(results_dir),
        subsample_n=50000,
        use_gpu="auto",
        max_epochs=200,
        save_intermediates=True,
    )

    logger.info("Integration benchmark results:")
    logger.info("\n%s", comparison_df.to_string())
    logger.info("Selected method: %s", best_method)

    # Save benchmark results
    comparison_df.to_csv(results_dir / "integration_benchmark.csv", index=False)

    # Transfer embeddings from HVG adata back to full adata
    for key in adata_hvg.obsm:
        if key not in ("spatial",):
            adata.obsm[key] = adata_hvg.obsm[key]
            logger.info("Transferred %s to full adata", key)

    # Neighbors, UMAP, Leiden on the best embedding
    neighbors(adata)  # auto-detects best use_rep
    umap(adata)
    cluster(adata, resolution=1.0)
    logger.info("Clustering complete: %d leiden clusters", adata.obs["leiden"].nunique())

    # Generate post-integration QC report
    figures_dir = PROJECT_ROOT / "figures"
    try:
        from sc_tools.qc.report import generate_post_integration_report

        generate_post_integration_report(
            adata,
            output_dir=str(figures_dir / "QC"),
            comparison_df=comparison_df,
            batch_key="sample",
        )
        logger.info("Post-integration QC report saved")
    except Exception as exc:
        logger.warning("Post-integration report failed (non-fatal): %s", exc)

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)
    logger.info(
        "Saved: %s (%d cells x %d genes, %d clusters)",
        output_path, adata.n_obs, adata.n_vars, adata.obs["leiden"].nunique(),
    )


if __name__ == "__main__":
    main()
