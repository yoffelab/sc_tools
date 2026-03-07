#!/usr/bin/env python3
"""
Figure 3: Clinical associations of LME classes.

Manuscript caption (adapted):
  a) KM curves for OS by LME class (n=268, log-rank p=0.0064)
  b) KM curves for PFS by LME class (n=264, log-rank p=0.0025)
  c) Cox forest plot (HR=2.69 Cytotoxic vs CD206 Enriched reference)
  d) COO x LME enrichment (Fisher exact + BH)
  e) Mutation frequency per LME class

Insight: LME classes have distinct clinical outcomes. Cytotoxic LME has the
worst overall and progression-free survival. CD206 Enriched is protective.
COO subtypes are non-uniformly distributed across LMEs.

Usage:
    python scripts/fig3_clinical.py

Input:
    metadata/DLC380_clinical.tsv
    metadata/lme_class_assignments.csv

Output:
    figures/manuscript/fig3/
"""

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

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_config import (
    COO_COLORS,
    LME_COLORS,
    LME_ORDER,
    apply_figure_style,
    significance_label,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

PROJECT_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_DIR / "results"
FIG_DIR = PROJECT_DIR / "figures" / "manuscript" / "fig3"
METADATA_DIR = PROJECT_DIR / "metadata"


def load_config():
    with open(PROJECT_DIR / "config.yaml") as f:
        return yaml.safe_load(f)


def load_clinical_with_lme() -> pd.DataFrame | None:
    """Load DLC380_clinical.tsv joined with LME assignments."""
    clinical_path = METADATA_DIR / "DLC380_clinical.tsv"
    lme_path = METADATA_DIR / "lme_class_assignments.csv"

    if not clinical_path.exists():
        logger.error(f"Clinical file not found: {clinical_path}")
        return None

    clinical = pd.read_csv(clinical_path, sep="\t")
    logger.info(f"Clinical: {clinical.shape[0]} cases")

    # Filter final cohort
    if "FINAL_COHORT" in clinical.columns:
        clinical = clinical[clinical["FINAL_COHORT"] == "YES"].copy()
        logger.info(f"  Final cohort: {len(clinical)} cases")

    # Standardize column names
    rename = {
        "Overall survival (y)": "OS_time",
        "CODE_OS": "OS_event",
        "Progression free survival (y)": "PFS_time",
        "CODE_PFS": "PFS_event",
        "LYMPH2CX_COO": "COO",
        "AGE": "age",
        "SEX": "sex",
        "IPI": "IPI",
        "STAGE": "stage",
    }
    clinical = clinical.rename(columns={k: v for k, v in rename.items() if k in clinical.columns})
    clinical["DLC_ID"] = clinical["DLC_ID"].astype(str)

    # Join with LME (normalize DLC ID format: DLC0002 <-> DLC_0002)
    if lme_path.exists():
        import re

        def norm_dlc(s):
            m = re.match(r"DLC[_\s-]?(\d+)", str(s), re.IGNORECASE)
            return f"DLC_{int(m.group(1)):04d}" if m else str(s)

        lme = pd.read_csv(lme_path)
        lme["sample_norm"] = lme["sample"].apply(norm_dlc)
        clinical["DLC_ID_norm"] = clinical["DLC_ID"].apply(norm_dlc)

        merged = clinical.merge(lme, left_on="DLC_ID_norm", right_on="sample_norm", how="inner")
        logger.info(f"  Joined clinical + LME: {len(merged)} cases")

        # Validate all 5 classes
        present = set(merged["LME_class"].unique())
        missing = set(LME_ORDER) - present
        if missing:
            logger.warning(f"  Missing LME classes: {missing}")

        return merged

    logger.warning("LME assignments not found")
    return clinical


def add_number_at_risk_table(ax, kmfs, time_points, lme_classes):
    """Add number-at-risk table below KM plot."""
    # Create sub-axes below the main plot
    table_ax = ax.get_figure().add_axes(
        [ax.get_position().x0, ax.get_position().y0 - 0.15,
         ax.get_position().width, 0.12]
    )
    table_ax.set_xlim(ax.get_xlim())
    table_ax.set_ylim(0, len(lme_classes))
    table_ax.set_axis_off()

    for i, (lme, kmf) in enumerate(zip(lme_classes, kmfs, strict=False)):
        color = LME_COLORS.get(lme, "#999999")
        table_ax.text(-0.02, len(lme_classes) - i - 0.5, lme,
                      ha="right", va="center", fontsize=5, color=color,
                      transform=table_ax.get_yaxis_transform())
        for t in time_points:
            n_at_risk = (kmf.survival_function_.index >= t).sum()
            table_ax.text(t, len(lme_classes) - i - 0.5, str(n_at_risk),
                          ha="center", va="center", fontsize=5)


def fig3a_km_os(df: pd.DataFrame):
    """KM curves for overall survival by LME class.

    Insight: Cytotoxic LME has worst OS; CD206 Enriched is protective.
    Direction: Cytotoxic curve should be BELOW others.
    """
    logger.info("  Fig 3a: KM OS")
    try:
        from lifelines import KaplanMeierFitter
        from lifelines.statistics import multivariate_logrank_test
    except ImportError:
        logger.warning("  lifelines not installed; skipping KM plots")
        return

    if "OS_time" not in df.columns or "OS_event" not in df.columns:
        logger.warning("  OS columns not found")
        return

    data = df[["OS_time", "OS_event", "LME_class"]].dropna()
    data["OS_time"] = pd.to_numeric(data["OS_time"], errors="coerce")
    data["OS_event"] = pd.to_numeric(data["OS_event"], errors="coerce")
    data = data.dropna()

    if len(data) < 10 or data["LME_class"].nunique() < 2:
        logger.warning(f"  Insufficient OS data: {len(data)} cases, {data['LME_class'].nunique()} classes")
        return

    fig, ax = plt.subplots(figsize=(6, 5))

    lme_present = [l for l in LME_ORDER if l in data["LME_class"].unique()]
    kmfs = []

    for lme in lme_present:
        mask = data["LME_class"] == lme
        n = mask.sum()
        kmf_i = KaplanMeierFitter()
        kmf_i.fit(data.loc[mask, "OS_time"], data.loc[mask, "OS_event"],
                  label=f"{lme} (n={n})")
        kmf_i.plot_survival_function(ax=ax, color=LME_COLORS.get(lme), linewidth=1.2)
        kmfs.append(kmf_i)

    # Global multivariate log-rank test
    result = multivariate_logrank_test(data["OS_time"], data["LME_class"], data["OS_event"])
    p_val = result.p_value

    ax.text(0.95, 0.95, f"log-rank p = {p_val:.4f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=8,
            bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.8})

    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Overall Survival Probability")
    ax.set_title(f"Overall Survival by LME Class (n={len(data)})")
    ax.legend(loc="lower left", fontsize=7)
    ax.set_ylim(0, 1.05)
    sns.despine()

    out = FIG_DIR / "fig3a_km_os.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out} (p={p_val:.4f})")


def fig3b_km_pfs(df: pd.DataFrame):
    """KM curves for PFS by LME class.

    Insight: Same pattern as OS — Cytotoxic worst, CD206 Enriched best.
    """
    logger.info("  Fig 3b: KM PFS")
    try:
        from lifelines import KaplanMeierFitter
        from lifelines.statistics import multivariate_logrank_test
    except ImportError:
        return

    if "PFS_time" not in df.columns or "PFS_event" not in df.columns:
        logger.warning("  PFS columns not found")
        return

    data = df[["PFS_time", "PFS_event", "LME_class"]].dropna()
    data["PFS_time"] = pd.to_numeric(data["PFS_time"], errors="coerce")
    data["PFS_event"] = pd.to_numeric(data["PFS_event"], errors="coerce")
    data = data.dropna()

    if len(data) < 10 or data["LME_class"].nunique() < 2:
        logger.warning(f"  Insufficient PFS data: {len(data)} cases")
        return

    fig, ax = plt.subplots(figsize=(6, 5))
    lme_present = [l for l in LME_ORDER if l in data["LME_class"].unique()]

    for lme in lme_present:
        mask = data["LME_class"] == lme
        n = mask.sum()
        kmf = KaplanMeierFitter()
        kmf.fit(data.loc[mask, "PFS_time"], data.loc[mask, "PFS_event"],
                label=f"{lme} (n={n})")
        kmf.plot_survival_function(ax=ax, color=LME_COLORS.get(lme), linewidth=1.2)

    result = multivariate_logrank_test(data["PFS_time"], data["LME_class"], data["PFS_event"])
    p_val = result.p_value

    ax.text(0.95, 0.95, f"log-rank p = {p_val:.4f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=8,
            bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.8})

    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Progression-Free Survival Probability")
    ax.set_title(f"PFS by LME Class (n={len(data)})")
    ax.legend(loc="lower left", fontsize=7)
    ax.set_ylim(0, 1.05)
    sns.despine()

    out = FIG_DIR / "fig3b_km_pfs.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out} (p={p_val:.4f})")


def fig3c_cox_forest(df: pd.DataFrame):
    """Cox PH forest plot — HR with 95% CI.

    Insight: Cytotoxic LME has HR=2.69 vs CD206 Enriched (reference).
    Direction: Cytotoxic bar should extend RIGHT of HR=1.
    """
    logger.info("  Fig 3c: Cox forest plot")
    try:
        from lifelines import CoxPHFitter
    except ImportError:
        logger.warning("  lifelines not installed")
        return

    if "OS_time" not in df.columns or "OS_event" not in df.columns:
        return

    data = df[["LME_class", "OS_time", "OS_event"]].dropna()
    data["OS_time"] = pd.to_numeric(data["OS_time"], errors="coerce")
    data["OS_event"] = pd.to_numeric(data["OS_event"], errors="coerce")
    data = data.dropna()

    # Dummy encode with CD206 Enriched as reference (protective)
    dummies = pd.get_dummies(data["LME_class"], prefix="LME", drop_first=False)
    ref_col = "LME_CD206 Enriched"
    if ref_col in dummies.columns:
        dummies = dummies.drop(columns=[ref_col])
    else:
        dummies = dummies.iloc[:, 1:]  # drop first alphabetically

    cox_data = pd.concat([dummies, data[["OS_time", "OS_event"]]], axis=1)
    cox_data.columns = cox_data.columns.str.replace(" ", "_")

    cph = CoxPHFitter()
    try:
        cph.fit(cox_data, duration_col="OS_time", event_col="OS_event")
        summary = cph.summary

        # Manual forest plot
        fig, ax = plt.subplots(figsize=(7, max(3, len(summary) * 0.5 + 1)))

        names = [n.replace("LME_", "").replace("_", " ") for n in summary.index]
        hrs = np.exp(summary["coef"])
        ci_low = np.exp(summary["coef"] - 1.96 * summary["se(coef)"])
        ci_high = np.exp(summary["coef"] + 1.96 * summary["se(coef)"])
        p_vals = summary["p"]

        y_pos = range(len(names))
        ax.errorbar(hrs, y_pos, xerr=[hrs - ci_low, ci_high - hrs],
                     fmt="o", color="black", capsize=3, markersize=5)
        ax.axvline(x=1, color="gray", linestyle="--", alpha=0.5, linewidth=0.8)
        ax.set_yticks(list(y_pos))
        ax.set_yticklabels(names)
        ax.set_xlabel("Hazard Ratio (95% CI)")
        ax.set_title("Cox PH — OS by LME Class\n(Reference: CD206 Enriched)")

        # Add HR and p-value labels
        for i, (hr, p) in enumerate(zip(hrs, p_vals, strict=False)):
            sig = significance_label(p)
            ax.text(max(ci_high) * 1.05, i,
                    f"HR={hr:.2f} (p={p:.3f}) {sig}",
                    va="center", fontsize=6)

        sns.despine()

        out = FIG_DIR / "fig3c_cox_forest.pdf"
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"    Saved: {out}")

        summary.to_csv(FIG_DIR / "fig3c_cox_summary.csv")
    except Exception as e:
        logger.warning(f"  Cox model failed: {e}")


def fig3d_coo(df: pd.DataFrame):
    """COO x LME enrichment (Fisher exact + BH).

    Insight: LME classes show non-random association with COO subtypes.
    """
    logger.info("  Fig 3d: COO distribution")

    if "COO" not in df.columns:
        logger.warning("  No COO column")
        return

    if "LME_class" not in df.columns:
        return

    data = df.dropna(subset=["COO", "LME_class"])
    ct = pd.crosstab(data["LME_class"], data["COO"])

    # Reorder
    lme_present = [l for l in LME_ORDER if l in ct.index]
    ct = ct.loc[lme_present]

    # Chi-squared test
    chi2, p_val, dof, expected = stats.chi2_contingency(ct)

    ct_norm = ct.div(ct.sum(axis=1), axis=0)
    coo_cols = ct_norm.columns.tolist()
    colors = [COO_COLORS.get(c, "#999999") for c in coo_cols]

    fig, ax = plt.subplots(figsize=(7, 5))
    ct_norm.plot(kind="bar", stacked=True, ax=ax, width=0.8, color=colors)
    ax.set_ylabel("Proportion")
    ax.set_xlabel("")
    ax.set_title(f"Cell of Origin by LME Class (chi-squared p = {p_val:.2e})")
    ax.legend(title="COO", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.xticks(rotation=30, ha="right")

    # Add n labels on bars
    for i, lme in enumerate(lme_present):
        n = ct.loc[lme].sum()
        ax.text(i, 1.02, f"n={n}", ha="center", fontsize=6)

    sns.despine()

    out = FIG_DIR / "fig3d_coo.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig3e_mutations(config: dict, clinical_df: pd.DataFrame):
    """Mutation frequency per LME class.

    Insight: Certain mutations are enriched in specific LME classes.
    """
    logger.info("  Fig 3e: Mutation frequency")

    mut_path = PROJECT_DIR / config["clinical"]["mutation_table"]
    if not mut_path.exists():
        logger.warning(f"  Mutation table not found: {mut_path}")
        return

    mut = pd.read_csv(mut_path)
    if "LME_class" not in clinical_df.columns:
        return

    # Join mutations with LME
    for col in mut.columns:
        overlap = set(mut[col].astype(str)) & set(clinical_df.get("DLC_ID", clinical_df.get("sample", pd.Series())).astype(str))
        if len(overlap) > 10:
            merged = mut.merge(
                clinical_df[["DLC_ID", "LME_class"]].drop_duplicates() if "DLC_ID" in clinical_df.columns
                else clinical_df[["sample", "LME_class"]].drop_duplicates(),
                left_on=col,
                right_on="DLC_ID" if "DLC_ID" in clinical_df.columns else "sample",
                how="inner",
            )
            break
    else:
        logger.warning("  Could not join mutations with LME")
        return

    gene_cols = [c for c in merged.columns
                 if c not in ["sample", "LME_class", col, "DLC_ID"]]
    if not gene_cols:
        return

    # Mean mutation rate per LME
    mut_rates = merged.groupby("LME_class")[gene_cols].mean()
    lme_present = [l for l in LME_ORDER if l in mut_rates.index]
    mut_rates = mut_rates.loc[lme_present]

    # Select genes with highest variance across LMEs
    gene_var = mut_rates.var(axis=0).sort_values(ascending=False)
    top_genes = gene_var.head(min(15, len(gene_var))).index.tolist()

    # Heatmap
    fig, ax = plt.subplots(figsize=(max(6, len(top_genes) * 0.4), 4))
    sns.heatmap(
        mut_rates[top_genes],
        cmap="YlOrRd",
        ax=ax,
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"label": "Mutation Rate"},
        linewidths=0.5,
    )
    ax.set_title("Mutation Frequency by LME Class")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    sns.despine()

    out = FIG_DIR / "fig3e_mutations.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def main():
    apply_figure_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    config = load_config()

    df = load_clinical_with_lme()
    if df is None:
        logger.error("Cannot generate Figure 3 without clinical data")
        return

    logger.info(f"Clinical columns: {df.columns.tolist()[:20]}")

    fig3a_km_os(df)
    fig3b_km_pfs(df)
    fig3c_cox_forest(df)
    fig3d_coo(df)
    fig3e_mutations(config, df)

    logger.info("Figure 3 complete.")


if __name__ == "__main__":
    main()
