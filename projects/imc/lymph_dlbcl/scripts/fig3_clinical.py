#!/usr/bin/env python3
"""
Figure 3: Clinical analysis — survival, COO, mutations.

Panels:
  a) KM curves for OS by LME class
  b) KM curves for PFS by LME class
  c) COO distribution per LME class (stacked barplot + chi-squared)
  d) Mutation frequency per LME class
  e) Forest plot of Cox regression hazard ratios

Usage:
    python scripts/fig3_clinical.py

Input:
    results/adata.*.celltyped.p4.h5ad (for LME + sample mapping)
    data/downloaded/clinical/DLBCL_clinical_full.csv
    data/downloaded/clinical/CTMA121_mut_table.csv
    metadata/lme_class_assignments.csv

Output:
    figures/manuscript/fig3/
"""

import logging
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
RESULTS_DIR = PROJECT_DIR / "results"
FIG_DIR = PROJECT_DIR / "figures" / "manuscript" / "fig3"
METADATA_DIR = PROJECT_DIR / "metadata"

LME_COLORS = {
    "Cold": "#4575b4",
    "CD206 Enriched": "#d73027",
    "Cytotoxic": "#fc8d59",
    "Stromal": "#91bfdb",
    "T cell Regulated": "#fee090",
}


def load_config():
    with open(PROJECT_DIR / "config.yaml") as f:
        return yaml.safe_load(f)


def load_clinical_with_lme(config: dict) -> pd.DataFrame | None:
    """Load clinical data joined with LME assignments."""
    clinical_path = PROJECT_DIR / config["clinical"]["clinical_full"]
    lme_path = METADATA_DIR / "lme_class_assignments.csv"

    if not clinical_path.exists():
        logger.warning(f"Clinical file not found: {clinical_path}")
        return None

    clinical = pd.read_csv(clinical_path)
    logger.info(f"Clinical: {clinical.shape}, columns: {clinical.columns.tolist()[:15]}")

    if lme_path.exists():
        lme = pd.read_csv(lme_path)
        lme["sample"] = lme["sample"].astype(str)

        # Find join key
        for col in clinical.columns:
            overlap = set(clinical[col].astype(str)) & set(lme["sample"])
            if len(overlap) > 10:
                clinical[col] = clinical[col].astype(str)
                merged = clinical.merge(
                    lme[["sample", "LME_display"]].rename(columns={"LME_display": "LME_class"}),
                    left_on=col, right_on="sample", how="inner",
                )
                logger.info(f"Joined clinical + LME via '{col}': {len(merged)} samples")
                return merged

    logger.warning("Could not join clinical with LME")
    return clinical


def find_survival_cols(df: pd.DataFrame) -> dict:
    """Auto-detect survival column names."""
    cols = {}
    for col in df.columns:
        cl = col.lower()
        if "os" in cl and ("time" in cl or "month" in cl or "year" in cl):
            cols["os_time"] = col
        elif "os" in cl and ("event" in cl or "status" in cl or "censor" in cl):
            cols["os_event"] = col
        elif "pfs" in cl and ("time" in cl or "month" in cl or "year" in cl):
            cols["pfs_time"] = col
        elif "pfs" in cl and ("event" in cl or "status" in cl or "censor" in cl):
            cols["pfs_event"] = col
        elif "coo" in cl or "cell_of_origin" in cl:
            cols["coo"] = col
    return cols


def fig3a_km_os(df: pd.DataFrame, surv_cols: dict):
    """KM curves for overall survival by LME class."""
    logger.info("  Fig 3a: KM OS")
    try:
        from lifelines import KaplanMeierFitter
        from lifelines.statistics import logrank_test
    except ImportError:
        logger.warning("  lifelines not installed; skipping KM plots")
        return

    if "os_time" not in surv_cols or "os_event" not in surv_cols:
        logger.warning("  OS columns not found")
        return

    time_col = surv_cols["os_time"]
    event_col = surv_cols["os_event"]

    data = df.dropna(subset=[time_col, event_col, "LME_class"])
    data[time_col] = pd.to_numeric(data[time_col], errors="coerce")
    data[event_col] = pd.to_numeric(data[event_col], errors="coerce")
    data = data.dropna(subset=[time_col, event_col])

    fig, ax = plt.subplots(figsize=(8, 6))
    kmf = KaplanMeierFitter()

    lme_classes = sorted(data["LME_class"].unique())
    for lme in lme_classes:
        mask = data["LME_class"] == lme
        kmf.fit(
            data.loc[mask, time_col],
            data.loc[mask, event_col],
            label=f"{lme} (n={mask.sum()})",
        )
        kmf.plot_survival_function(ax=ax, color=LME_COLORS.get(lme, None))

    # Global log-rank test (pairwise)
    if len(lme_classes) >= 2:
        p_values = []
        for i, lme1 in enumerate(lme_classes):
            for lme2 in lme_classes[i + 1:]:
                m1 = data["LME_class"] == lme1
                m2 = data["LME_class"] == lme2
                result = logrank_test(
                    data.loc[m1, time_col], data.loc[m2, time_col],
                    data.loc[m1, event_col], data.loc[m2, event_col],
                )
                p_values.append(result.p_value)
        if p_values:
            min_p = min(p_values)
            ax.text(0.95, 0.95, f"min log-rank p = {min_p:.2e}",
                    transform=ax.transAxes, ha="right", va="top", fontsize=9)

    ax.set_xlabel("Time")
    ax.set_ylabel("Overall Survival Probability")
    ax.set_title("Overall Survival by LME Class")
    ax.legend(loc="lower left", fontsize=8)
    plt.tight_layout()

    out = FIG_DIR / "fig3a_km_os.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig3b_km_pfs(df: pd.DataFrame, surv_cols: dict):
    """KM curves for PFS by LME class."""
    logger.info("  Fig 3b: KM PFS")
    try:
        from lifelines import KaplanMeierFitter
    except ImportError:
        return

    if "pfs_time" not in surv_cols or "pfs_event" not in surv_cols:
        logger.warning("  PFS columns not found")
        return

    time_col = surv_cols["pfs_time"]
    event_col = surv_cols["pfs_event"]

    data = df.dropna(subset=[time_col, event_col, "LME_class"])
    data[time_col] = pd.to_numeric(data[time_col], errors="coerce")
    data[event_col] = pd.to_numeric(data[event_col], errors="coerce")
    data = data.dropna(subset=[time_col, event_col])

    fig, ax = plt.subplots(figsize=(8, 6))
    kmf = KaplanMeierFitter()

    for lme in sorted(data["LME_class"].unique()):
        mask = data["LME_class"] == lme
        kmf.fit(data.loc[mask, time_col], data.loc[mask, event_col],
                label=f"{lme} (n={mask.sum()})")
        kmf.plot_survival_function(ax=ax, color=LME_COLORS.get(lme, None))

    ax.set_xlabel("Time")
    ax.set_ylabel("Progression-Free Survival Probability")
    ax.set_title("PFS by LME Class")
    ax.legend(loc="lower left", fontsize=8)
    plt.tight_layout()

    out = FIG_DIR / "fig3b_km_pfs.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig3c_coo(df: pd.DataFrame, surv_cols: dict):
    """COO distribution per LME class."""
    logger.info("  Fig 3c: COO distribution")

    if "coo" not in surv_cols:
        # Try to find COO column
        coo_col = None
        for col in df.columns:
            if "coo" in col.lower() or "abc" in str(df[col].unique()).lower():
                coo_col = col
                break
        if coo_col is None:
            logger.warning("  No COO column found")
            return
    else:
        coo_col = surv_cols["coo"]

    if "LME_class" not in df.columns:
        return

    data = df.dropna(subset=[coo_col, "LME_class"])
    ct = pd.crosstab(data["LME_class"], data[coo_col])

    # Chi-squared test
    chi2, p_val, dof, expected = stats.chi2_contingency(ct)

    fig, ax = plt.subplots(figsize=(10, 6))
    ct_norm = ct.div(ct.sum(axis=1), axis=0)
    ct_norm.plot(kind="bar", stacked=True, ax=ax, width=0.8)
    ax.set_ylabel("Proportion")
    ax.set_xlabel("")
    ax.set_title(f"Cell of Origin by LME Class (chi2 p = {p_val:.2e})")
    ax.legend(title="COO", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()

    out = FIG_DIR / "fig3c_coo.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig3d_mutations(config: dict, clinical_df: pd.DataFrame):
    """Mutation frequency per LME class."""
    logger.info("  Fig 3d: Mutation frequency")

    mut_path = PROJECT_DIR / config["clinical"]["mutation_table"]
    if not mut_path.exists():
        logger.warning(f"  Mutation table not found: {mut_path}")
        return

    mut = pd.read_csv(mut_path)
    logger.info(f"  Mutation table: {mut.shape}")

    if "LME_class" not in clinical_df.columns:
        return

    # Join mutations with LME via sample ID
    for col in mut.columns:
        overlap = set(mut[col].astype(str)) & set(clinical_df.get("sample", pd.Series()).astype(str))
        if len(overlap) > 10:
            merged = mut.merge(
                clinical_df[["sample", "LME_class"]].drop_duplicates(),
                left_on=col, right_on="sample", how="inner",
            )
            break
    else:
        logger.warning("  Could not join mutations with LME")
        return

    # Compute mutation rates per LME
    gene_cols = [c for c in merged.columns if c not in ["sample", "LME_class", col]]
    if not gene_cols:
        return

    mut_rates = merged.groupby("LME_class")[gene_cols].mean()

    fig, ax = plt.subplots(figsize=(max(10, len(gene_cols) * 0.5), 6))
    mut_rates.T.plot(kind="bar", ax=ax, width=0.8)
    ax.set_ylabel("Mutation Rate")
    ax.set_xlabel("Gene")
    ax.set_title("Mutation Frequency by LME Class")
    ax.legend(title="LME", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    plt.xticks(rotation=90, fontsize=8)
    plt.tight_layout()

    out = FIG_DIR / "fig3d_mutations.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig3e_cox_forest(df: pd.DataFrame, surv_cols: dict):
    """Forest plot of Cox regression hazard ratios."""
    logger.info("  Fig 3e: Cox forest plot")
    try:
        from lifelines import CoxPHFitter
    except ImportError:
        logger.warning("  lifelines not installed")
        return

    if "os_time" not in surv_cols or "os_event" not in surv_cols:
        return

    time_col = surv_cols["os_time"]
    event_col = surv_cols["os_event"]

    data = df[["LME_class", time_col, event_col]].dropna()
    data[time_col] = pd.to_numeric(data[time_col], errors="coerce")
    data[event_col] = pd.to_numeric(data[event_col], errors="coerce")
    data = data.dropna()

    # Dummy encode LME (reference = Cold)
    dummies = pd.get_dummies(data["LME_class"], prefix="LME", drop_first=True)
    cox_data = pd.concat([dummies, data[[time_col, event_col]]], axis=1)
    cox_data.columns = cox_data.columns.str.replace(" ", "_")

    cph = CoxPHFitter()
    try:
        cph.fit(cox_data, duration_col=time_col, event_col=event_col)

        fig, ax = plt.subplots(figsize=(8, 5))
        cph.plot(ax=ax)
        ax.set_title("Cox Proportional Hazards — OS by LME Class\n(Reference: Cold)")
        ax.axvline(x=0, color="gray", linestyle="--", alpha=0.5)
        plt.tight_layout()

        out = FIG_DIR / "fig3e_cox_forest.pdf"
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"    Saved: {out}")

        # Save summary
        summary = cph.summary
        summary.to_csv(FIG_DIR / "fig3e_cox_summary.csv")
    except Exception as e:
        logger.warning(f"  Cox model failed: {e}")


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    config = load_config()

    df = load_clinical_with_lme(config)
    if df is None:
        logger.error("Cannot generate Figure 3 without clinical data")
        return

    surv_cols = find_survival_cols(df)
    logger.info(f"Detected survival columns: {surv_cols}")

    fig3a_km_os(df, surv_cols)
    fig3b_km_pfs(df, surv_cols)
    fig3c_coo(df, surv_cols)
    fig3d_mutations(config, df)
    fig3e_cox_forest(df, surv_cols)

    logger.info("Figure 3 complete.")


if __name__ == "__main__":
    main()
