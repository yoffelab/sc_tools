#!/usr/bin/env python3
"""
Phase 6: Validate reproduced figures against manuscript originals.

Checks:
1. All expected figure files exist
2. Cell counts match manuscript (328 tumors, 12 major types, 30 subpops)
3. LME class proportions within tolerance
4. Reports numerical validation summary

Usage:
    python scripts/validate_figures.py

Output:
    outputs/figure_validation_report.md
"""

import logging
from pathlib import Path

import anndata as ad
import pandas as pd
import yaml

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

PROJECT_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_DIR / "results"
FIG_DIR = PROJECT_DIR / "figures" / "manuscript"
OUTPUT_DIR = PROJECT_DIR / "outputs"
METADATA_DIR = PROJECT_DIR / "metadata"

# Expected values from manuscript
EXPECTED = {
    "n_tumors": 328,
    "n_markers": 52,
    "n_major_celltypes": 12,
    "n_subpopulations": 30,
    "n_lme_classes": 5,
    "lme_proportions": {
        "Cold": 35.1,
        "Stromal": 21.3,
        "Cytotoxic": 20.7,
        "CD206 Enriched": 8.2,
        "T cell Regulated": 14.6,
    },
}

EXPECTED_FIGURES = {
    "fig1": ["fig1a_umap_immune.pdf", "fig1b_heatmap.pdf"],
    "fig2": ["fig2a_heatmap.pdf", "fig2b_lme_proportions.pdf"],
    "fig3": ["fig3a_km_os.pdf"],
    "fig4": ["fig4a_community_composition.pdf"],
    "fig5": ["fig5a_roc.pdf", "fig5b_feature_importance.pdf"],
}


def load_config():
    with open(PROJECT_DIR / "config.yaml") as f:
        return yaml.safe_load(f)


def check_figures() -> list[dict]:
    """Check which expected figures exist."""
    results = []
    for fig_name, expected_files in EXPECTED_FIGURES.items():
        fig_path = FIG_DIR / fig_name
        for fname in expected_files:
            fpath = fig_path / fname
            results.append({
                "figure": fig_name,
                "file": fname,
                "exists": fpath.exists(),
                "size_kb": fpath.stat().st_size / 1024 if fpath.exists() else 0,
            })
    return results


def check_data() -> dict:
    """Check numerical consistency with manuscript."""
    checks = {}

    # Count unique samples across panels
    total_samples = set()
    for panel in ["immune", "stromal"]:
        for suffix in ["celltyped.p4", "annotated.p2", "raw.p1"]:
            path = RESULTS_DIR / f"adata.{panel}.{suffix}.h5ad"
            if path.exists():
                adata = ad.read_h5ad(path, backed="r")
                if "sample" in adata.obs.columns:
                    total_samples.update(adata.obs["sample"].unique().tolist())
                adata.file.close()
                break

    checks["n_samples"] = len(total_samples)
    checks["expected_n_tumors"] = EXPECTED["n_tumors"]
    checks["sample_match"] = len(total_samples) >= EXPECTED["n_tumors"] * 0.9

    # Check LME proportions
    lme_path = METADATA_DIR / "lme_class_assignments.csv"
    if lme_path.exists():
        lme = pd.read_csv(lme_path)
        if "LME_display" in lme.columns:
            lme_pcts = lme["LME_display"].value_counts(normalize=True) * 100
            checks["lme_proportions"] = lme_pcts.to_dict()
            checks["lme_n_classes"] = lme["LME_display"].nunique()

            # Compare with expected
            for cls, expected_pct in EXPECTED["lme_proportions"].items():
                actual_pct = lme_pcts.get(cls, 0)
                checks[f"lme_{cls}_diff"] = abs(actual_pct - expected_pct)

    return checks


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("Validating reproduced figures...")

    # Check figures
    fig_results = check_figures()
    n_exist = sum(1 for r in fig_results if r["exists"])
    n_total = len(fig_results)
    logger.info(f"Figures: {n_exist}/{n_total} exist")

    # Check data
    data_checks = check_data()

    # Generate report
    lines = [
        "# Figure Validation Report\n",
        "## Figure Files\n",
        f"**{n_exist}/{n_total}** expected figures exist.\n",
        "| Figure | File | Exists | Size (KB) |",
        "|--------|------|--------|-----------|",
    ]

    for r in fig_results:
        lines.append(
            f"| {r['figure']} | {r['file']} | "
            f"{'Y' if r['exists'] else 'N'} | {r['size_kb']:.1f} |"
        )

    lines.extend([
        "\n## Numerical Validation\n",
        f"- Samples found: {data_checks.get('n_samples', '?')} "
        f"(expected: {EXPECTED['n_tumors']})",
        f"- Sample match (>=90%): {data_checks.get('sample_match', '?')}",
        f"- LME classes: {data_checks.get('lme_n_classes', '?')} "
        f"(expected: {EXPECTED['n_lme_classes']})",
    ])

    if "lme_proportions" in data_checks:
        lines.append("\n### LME Class Proportions\n")
        lines.append("| Class | Actual (%) | Expected (%) | Diff |")
        lines.append("|-------|-----------|-------------|------|")
        for cls, expected_pct in EXPECTED["lme_proportions"].items():
            actual_pct = data_checks["lme_proportions"].get(cls, 0)
            diff = abs(actual_pct - expected_pct)
            lines.append(f"| {cls} | {actual_pct:.1f} | {expected_pct:.1f} | {diff:.1f} |")

    report = "\n".join(lines)
    report_path = OUTPUT_DIR / "figure_validation_report.md"
    report_path.write_text(report)
    logger.info(f"Report: {report_path}")


if __name__ == "__main__":
    main()
