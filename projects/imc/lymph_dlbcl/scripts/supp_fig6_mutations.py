#!/usr/bin/env python3
"""Supp Fig 6: Mutation landscape — per-gene rates, co-occurrence, LME enrichment."""

import logging
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from scipy import stats
from statsmodels.stats.multitest import multipletests

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

PROJECT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_config import apply_figure_style

FIG_DIR = PROJECT_DIR / "figures" / "manuscript" / "supp_fig6"
METADATA_DIR = PROJECT_DIR / "metadata"


def load_config():
    with open(PROJECT_DIR / "config.yaml") as f:
        return yaml.safe_load(f)


def main():
    apply_figure_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    config = load_config()

    mut_path = PROJECT_DIR / config["clinical"]["mutation_table"]
    if not mut_path.exists():
        logger.error(f"Mutation table not found: {mut_path}")
        return

    mut = pd.read_csv(mut_path)
    logger.info(f"Mutation table: {mut.shape}")

    # Identify sample and gene columns
    sample_col = mut.columns[0]
    gene_cols = [c for c in mut.columns[1:] if mut[c].dtype in [int, float, np.int64, np.float64, bool]
                 or set(mut[c].dropna().unique()) <= {0, 1, "0", "1"}]

    if not gene_cols:
        gene_cols = mut.columns[1:].tolist()
        for col in gene_cols:
            mut[col] = pd.to_numeric(mut[col], errors="coerce").fillna(0).astype(int)

    logger.info(f"  {len(gene_cols)} genes, {len(mut)} samples")

    # a) Overall mutation rates
    mut_rates = mut[gene_cols].mean().sort_values(ascending=False)

    fig, ax = plt.subplots(figsize=(max(8, len(gene_cols) * 0.4), 5))
    ax.bar(range(len(mut_rates)), mut_rates.values, color="#4575b4")
    ax.set_xticks(range(len(mut_rates)))
    ax.set_xticklabels(mut_rates.index, rotation=90, fontsize=8)
    ax.set_ylabel("Mutation Rate")
    ax.set_title("Per-Gene Mutation Frequency")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "supp6a_mutation_rates.pdf", dpi=300, bbox_inches="tight")
    plt.close()

    # b) Co-occurrence matrix
    if len(gene_cols) > 1:
        cooccurrence = pd.DataFrame(index=gene_cols, columns=gene_cols, dtype=float)
        for g1 in gene_cols:
            for g2 in gene_cols:
                both = ((mut[g1] > 0) & (mut[g2] > 0)).sum()
                either = ((mut[g1] > 0) | (mut[g2] > 0)).sum()
                cooccurrence.loc[g1, g2] = both / max(either, 1)

        fig, ax = plt.subplots(figsize=(max(8, len(gene_cols) * 0.4), max(6, len(gene_cols) * 0.35)))
        sns.heatmap(cooccurrence.astype(float), cmap="YlOrRd", ax=ax,
                    xticklabels=True, yticklabels=True, vmin=0, vmax=1)
        ax.set_title("Mutation Co-occurrence (Jaccard)")
        plt.setp(ax.get_xticklabels(), fontsize=7, rotation=90)
        plt.setp(ax.get_yticklabels(), fontsize=7)
        plt.tight_layout()
        plt.savefig(FIG_DIR / "supp6b_cooccurrence.pdf", dpi=300, bbox_inches="tight")
        plt.close()

    # c) Mutation enrichment per LME class
    lme_path = METADATA_DIR / "lme_class_assignments.csv"
    if lme_path.exists():
        lme = pd.read_csv(lme_path)
        lme["sample"] = lme["sample"].astype(str)
        mut[sample_col] = mut[sample_col].astype(str)

        merged = mut.merge(lme[["sample", "LME_display"]], left_on=sample_col, right_on="sample", how="inner")
        if not merged.empty:
            lme_rates = merged.groupby("LME_display")[gene_cols].mean()

            fig, ax = plt.subplots(figsize=(max(10, len(gene_cols) * 0.4), 6))
            sns.heatmap(lme_rates, cmap="YlOrRd", ax=ax, vmin=0,
                        xticklabels=True, yticklabels=True,
                        cbar_kws={"label": "Mutation Rate"})
            ax.set_title("Mutation Rate per LME Class")
            plt.setp(ax.get_xticklabels(), fontsize=7, rotation=90)
            plt.tight_layout()
            plt.savefig(FIG_DIR / "supp6c_mutation_lme.pdf", dpi=300, bbox_inches="tight")
            plt.close()

            # Fisher exact test per gene x LME
            results = []
            for gene in gene_cols:
                for lme_cls in merged["LME_display"].unique():
                    in_lme = merged["LME_display"] == lme_cls
                    table = pd.crosstab(in_lme, merged[gene] > 0)
                    if table.shape == (2, 2):
                        _, p = stats.fisher_exact(table)
                        results.append({"gene": gene, "LME": lme_cls, "p_value": p})

            if results:
                res_df = pd.DataFrame(results)
                _, res_df["p_adj"], _, _ = multipletests(res_df["p_value"], method="fdr_bh")
                sig = res_df[res_df["p_adj"] < 0.05]
                if not sig.empty:
                    logger.info("  Significant mutation-LME associations (BH<0.05):")
                    for _, row in sig.iterrows():
                        logger.info(f"    {row['gene']} ~ {row['LME']}: p_adj={row['p_adj']:.3e}")
                res_df.to_csv(FIG_DIR / "supp6c_fisher_tests.csv", index=False)

    logger.info("Supp Fig 6 complete.")


if __name__ == "__main__":
    main()
