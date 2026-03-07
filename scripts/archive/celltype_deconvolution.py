"""
Cell-type deconvolution for GGO Visium project.

Thin wrapper around ``sc_tools.tl.deconvolution()``.
Run from project root:  python scripts/celltype_deconvolution.py [--method cell2location|tangram|destvi]

Default method: cell2location (requires GPU).
Use --method tangram for CPU/macOS fallback.
"""

import argparse
import logging
import sys

from sc_tools.tl.deconvolution import deconvolution

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cell-type deconvolution for GGO Visium")
    parser.add_argument("--method", default="cell2location", choices=["cell2location", "tangram", "destvi"])
    args = parser.parse_args()

    deconvolution(
        spatial_adata=(
            "results/adata.scored.h5ad"
            if __import__("os").path.exists("results/adata.scored.h5ad")
            else "results/adata.normalized.scored.p35.h5ad"
        ),
        sc_adata="results/seurat_object.h5ad",
        method=args.method,
        celltype_key="cell.type",
        spatial_batch_key="library_id",
        sc_batch_key="Batch",
        qc_labels=["QC_Filtered", "Doublets", "Low quality", "Unknown III (SM)"],
        n_signature_genes=2000,
        output_file=f"results/adata.deconvolution.{args.method}.h5ad",
        output_dir="output/deconvolution",
        cache_dir="output/deconvolution/cache",
    )
