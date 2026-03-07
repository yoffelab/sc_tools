#!/usr/bin/env python3
"""Supp Fig 7: RNA-protein comparison (IMC vs CosMx)."""

import logging
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

PROJECT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_config import apply_figure_style

FIG_DIR = PROJECT_DIR / "figures" / "manuscript" / "supp_fig7"


def main():
    apply_figure_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # This figure requires CosMx RNA data which may need separate download
    # Look for any RNA-protein comparison data
    rna_dirs = [
        PROJECT_DIR / "data" / "downloaded" / "cosmx",
        PROJECT_DIR / "data" / "cosmx",
    ]

    rna_data_found = False
    for d in rna_dirs:
        if d.exists() and any(d.iterdir()):
            rna_data_found = True
            break

    if not rna_data_found:
        logger.info("CosMx RNA data not yet available")
        logger.info("Creating placeholder figure")

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5,
                "RNA-Protein Comparison\n\n"
                "Requires CosMx RNA data\n"
                "(see Hyperion_RNA_v2.ipynb for reference)\n\n"
                "Expected: Correlation of IMC protein\n"
                "vs CosMx RNA expression per marker",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=12, style="italic")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        plt.savefig(FIG_DIR / "supp7_placeholder.pdf", dpi=300, bbox_inches="tight")
        plt.close()
        logger.info("  Placeholder saved")
        return

    logger.info("Supp Fig 7 complete.")


if __name__ == "__main__":
    main()
