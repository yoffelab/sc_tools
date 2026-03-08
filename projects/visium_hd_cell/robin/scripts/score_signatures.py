"""Phase 3.5b (scoring): Gene signature scoring for robin cell data.

Reads results/adata.normalized.h5ad and metadata/gene_signatures.json,
scores all signatures (project + Hallmark), saves results/adata.scored.h5ad.

Run from project root or via SLURM.
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

    from sc_tools.tl import score_signature
    from sc_tools.utils.checkpoint import resolve_checkpoint_path

    p3_path = resolve_checkpoint_path("preprocess", PROJECT_ROOT)
    if not p3_path.exists():
        logger.error("Preprocessed adata not found: %s", p3_path)
        sys.exit(1)

    sigs_path = PROJECT_ROOT / "metadata" / "gene_signatures.json"
    if not sigs_path.exists():
        logger.error("Gene signatures not found: %s", sigs_path)
        sys.exit(1)

    out_path = PROJECT_ROOT / "results" / "adata.scored.h5ad"

    logger.info("Reading %s", p3_path)
    adata = sc.read_h5ad(p3_path)
    adata.var_names_make_unique()
    logger.info("Loaded: %d cells x %d genes", adata.n_obs, adata.n_vars)

    logger.info("Scoring signatures from %s", sigs_path)
    score_signature(
        adata,
        signatures_nested=str(sigs_path),
        use_raw=True,
        ctrl_size=50,
        n_bins=25,
        min_genes=3,
        copy=False,
        save_path=str(out_path),
    )
    logger.info("Saved scored adata: %s", out_path)


if __name__ == "__main__":
    main()
