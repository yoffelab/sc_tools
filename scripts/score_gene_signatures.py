"""
Gene signature scoring for ggo_visium: load p2 adata and metadata/gene_signatures.json,
call sc_tools.tl.score_signature, write results/adata.normalized.scored.p35.h5ad.

Scores are stored in obsm (signature_score, signature_score_z) for traceability.
Run from project root: python scripts/score_gene_signatures.py
Requires sc_tools on PYTHONPATH (e.g. pip install -e . from repo root) or path fix below.

When running in Docker, use the project script by full path so the correct project dir
is used: python projects/visium/ggo_visium/scripts/score_gene_signatures.py
Alternatively set GGO_VISIUM_PROJECT_DIR to the project path (e.g. /workspace/projects/visium/ggo_visium).
"""

import os
import sys
from pathlib import Path

# Project root: env var for Docker/explicit override, else derive from this script's location
_GGO_ROOT_ENV = os.environ.get("GGO_VISIUM_PROJECT_DIR")
if _GGO_ROOT_ENV and Path(_GGO_ROOT_ENV).is_dir():
    GGO_ROOT = Path(_GGO_ROOT_ENV).resolve()
else:
    GGO_ROOT = Path(__file__).resolve().parent.parent

REPO_ROOT = GGO_ROOT.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import scanpy as sc

from sc_tools.tl import score_signature
ADATA_PATH = GGO_ROOT / "results" / "adata.annotated.p2.h5ad"
SIGNATURES_PATH = GGO_ROOT / "metadata" / "gene_signatures.json"
P35_PATH = GGO_ROOT / "results" / "adata.normalized.scored.p35.h5ad"


def main():
    if not ADATA_PATH.exists():
        raise FileNotFoundError(f"AnnData not found: {ADATA_PATH}. Run phase1 first.")
    if not SIGNATURES_PATH.exists():
        raise FileNotFoundError(f"Signatures not found: {SIGNATURES_PATH}")

    adata = sc.read_h5ad(ADATA_PATH)
    adata.var_names_make_unique()

    score_signature(
        adata,
        signatures_nested=str(SIGNATURES_PATH),
        use_raw=True,
        ctrl_size=50,
        n_bins=25,
        min_genes=3,
        copy=False,
        save_path=str(P35_PATH),
    )
    print(f"Saved: {P35_PATH}")


if __name__ == "__main__":
    main()
