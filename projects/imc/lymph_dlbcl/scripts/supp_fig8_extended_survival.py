#!/usr/bin/env python3
"""Supp Fig 8: Extended survival analysis and IHC validation."""

import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

PROJECT_DIR = Path(__file__).resolve().parent.parent
FIG_DIR = PROJECT_DIR / "figures" / "manuscript" / "supp_fig8"
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


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    config = load_config()

    try:
        from lifelines import CoxPHFitter, KaplanMeierFitter
        from lifelines.statistics import logrank_test
    except ImportError:
        logger.error("lifelines not installed; pip install lifelines")
        return

    # Load clinical + LME
    clinical_path = PROJECT_DIR / config["clinical"]["clinical_full"]
    lme_path = METADATA_DIR / "lme_class_assignments.csv"

    if not clinical_path.exists():
        logger.error(f"Clinical file not found: {clinical_path}")
        return

    clinical = pd.read_csv(clinical_path)

    if lme_path.exists():
        lme = pd.read_csv(lme_path)
        lme["sample"] = lme["sample"].astype(str)

        # Join
        for col in clinical.columns:
            overlap = set(clinical[col].astype(str)) & set(lme["sample"])
            if len(overlap) > 10:
                clinical[col] = clinical[col].astype(str)
                clinical = clinical.merge(
                    lme[["sample", "LME_display"]].rename(columns={"LME_display": "LME_class"}),
                    left_on=col, right_on="sample", how="inner",
                )
                break

    # Find survival columns
    surv_cols = {}
    for col in clinical.columns:
        cl = col.lower()
        if "os" in cl and ("time" in cl or "month" in cl):
            surv_cols["os_time"] = col
        elif "os" in cl and ("event" in cl or "status" in cl):
            surv_cols["os_event"] = col
        elif "pfs" in cl and ("time" in cl or "month" in cl):
            surv_cols["pfs_time"] = col
        elif "pfs" in cl and ("event" in cl or "status" in cl):
            surv_cols["pfs_event"] = col

    if "os_time" not in surv_cols:
        logger.error("No OS time column found")
        return

    time_col = surv_cols["os_time"]
    event_col = surv_cols.get("os_event", None)

    clinical[time_col] = pd.to_numeric(clinical[time_col], errors="coerce")
    if event_col:
        clinical[event_col] = pd.to_numeric(clinical[event_col], errors="coerce")
    clinical = clinical.dropna(subset=[time_col] + ([event_col] if event_col else []))

    # a) Multivariate Cox model
    if "LME_class" in clinical.columns:
        covariates = []
        for pattern in ["age", "sex", "stage", "ipi", "ecog", "ldh", "coo"]:
            for col in clinical.columns:
                if pattern in col.lower() and col not in covariates:
                    covariates.append(col)
                    break

        # LME dummies
        dummies = pd.get_dummies(clinical["LME_class"], prefix="LME", drop_first=True)

        # Numeric covariates only
        cox_cols = []
        for col in covariates:
            clinical[f"cov_{col}"] = pd.to_numeric(clinical[col], errors="coerce")
            if clinical[f"cov_{col}"].notna().sum() > len(clinical) * 0.5:
                cox_cols.append(f"cov_{col}")

        cox_data = pd.concat([
            dummies,
            clinical[cox_cols + [time_col] + ([event_col] if event_col else [])],
        ], axis=1).dropna()

        cox_data.columns = cox_data.columns.str.replace(" ", "_").str.replace(".", "_")
        time_col_clean = time_col.replace(" ", "_").replace(".", "_")
        event_col_clean = event_col.replace(" ", "_").replace(".", "_") if event_col else None

        if len(cox_data) > 20:
            cph = CoxPHFitter()
            try:
                cph.fit(
                    cox_data,
                    duration_col=time_col_clean,
                    event_col=event_col_clean,
                )

                fig, ax = plt.subplots(figsize=(10, max(5, len(cox_data.columns) * 0.3)))
                cph.plot(ax=ax)
                ax.axvline(x=0, color="gray", linestyle="--", alpha=0.5)
                ax.set_title("Multivariate Cox Regression — OS")
                plt.tight_layout()
                plt.savefig(FIG_DIR / "supp8a_multivariate_cox.pdf", dpi=300, bbox_inches="tight")
                plt.close()

                cph.summary.to_csv(FIG_DIR / "supp8a_cox_summary.csv")
                logger.info("  Multivariate Cox done")
            except Exception as e:
                logger.warning(f"  Cox model failed: {e}")

    # b) Subgroup KM curves (by COO)
    if "LME_class" in clinical.columns and event_col:
        coo_col = None
        for col in clinical.columns:
            if "coo" in col.lower():
                coo_col = col
                break

        if coo_col:
            coo_values = clinical[coo_col].dropna().unique()
            for coo in coo_values[:3]:  # Top 3 COO subtypes
                sub = clinical[clinical[coo_col] == coo].dropna(subset=[time_col, event_col])
                if len(sub) < 10:
                    continue

                fig, ax = plt.subplots(figsize=(8, 6))
                kmf = KaplanMeierFitter()
                for lme_cls in sorted(sub["LME_class"].unique()):
                    mask = sub["LME_class"] == lme_cls
                    if mask.sum() < 3:
                        continue
                    kmf.fit(sub.loc[mask, time_col], sub.loc[mask, event_col],
                            label=f"{lme_cls} (n={mask.sum()})")
                    kmf.plot_survival_function(ax=ax, color=LME_COLORS.get(lme_cls))

                ax.set_title(f"OS by LME — {coo} subgroup")
                ax.set_xlabel("Time")
                ax.set_ylabel("Survival Probability")
                ax.legend(fontsize=8)
                plt.tight_layout()

                coo_clean = str(coo).replace("/", "_").replace(" ", "_")
                plt.savefig(FIG_DIR / f"supp8b_km_{coo_clean}.pdf", dpi=300, bbox_inches="tight")
                plt.close()

    logger.info("Supp Fig 8 complete.")


if __name__ == "__main__":
    main()
