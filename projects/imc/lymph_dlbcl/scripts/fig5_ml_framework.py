#!/usr/bin/env python3
"""
Figure 5: ML framework for TME classification and IHC validation.

Panels:
  a) ROC curves for multiclass LME classification (one-vs-rest)
  b) Feature importance (top 20, horizontal barplot)
  c) Confusion matrix
  d) IHC validation scatter plot
  e) Performance metrics summary

Usage:
    python scripts/fig5_ml_framework.py

Input:
    data/downloaded/metadata/7.14.22.TME.zscore.csv
    data/downloaded/metadata/2.23.23.IHC_337_cases.csv
    data/downloaded/metadata/2.23.23.IHC_266_all_markers.csv
    metadata/lme_class_assignments.csv

Output:
    figures/manuscript/fig5/
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
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    classification_report,
    confusion_matrix,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import RandomizedSearchCV, StratifiedKFold
from sklearn.preprocessing import LabelBinarizer, StandardScaler

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

PROJECT_DIR = Path(__file__).resolve().parent.parent
FIG_DIR = PROJECT_DIR / "figures" / "manuscript" / "fig5"
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


def load_data(config: dict) -> tuple[pd.DataFrame, pd.Series] | tuple[None, None]:
    """Load TME z-scores + LME labels."""
    zscore_path = PROJECT_DIR / config["metadata"].get("tme_zscore", "")
    lme_path = METADATA_DIR / "lme_class_assignments.csv"

    if not zscore_path.exists():
        logger.warning(f"TME z-score file not found: {zscore_path}")
        return None, None
    if not lme_path.exists():
        logger.warning(f"LME assignments not found: {lme_path}")
        return None, None

    zscore = pd.read_csv(zscore_path, index_col=0)
    lme = pd.read_csv(lme_path)
    lme["sample"] = lme["sample"].astype(str)
    zscore.index = zscore.index.astype(str)

    merged = zscore.join(lme.set_index("sample")[["LME_display"]], how="inner")
    merged = merged.dropna(subset=["LME_display"])

    X = merged.drop(columns=["LME_display"])
    y = merged["LME_display"]

    logger.info(f"ML data: {X.shape[0]} samples, {X.shape[1]} features, {y.nunique()} classes")
    return X, y


def train_and_evaluate(X: pd.DataFrame, y: pd.Series, config: dict) -> dict:
    """Train RF with CV, return results."""
    seed = config["ml"]["random_seed"]
    n_splits = config["ml"]["n_splits"]
    n_iter = config["ml"]["n_iter"]

    scaler = StandardScaler()
    X_scaled = pd.DataFrame(scaler.fit_transform(X), index=X.index, columns=X.columns)

    # Hyperparameter search space
    param_dist = {
        "n_estimators": [100, 200, 500, 1000],
        "max_depth": [5, 10, 15, 20, None],
        "min_samples_split": [2, 5, 10],
        "min_samples_leaf": [1, 2, 4],
        "max_features": ["sqrt", "log2", None],
    }

    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    rf = RandomForestClassifier(random_state=seed, n_jobs=-1)

    search = RandomizedSearchCV(
        rf, param_dist, n_iter=n_iter, cv=cv, scoring="accuracy",
        random_state=seed, n_jobs=-1, return_train_score=True,
    )
    search.fit(X_scaled, y)

    best_model = search.best_estimator_
    logger.info(f"Best params: {search.best_params_}")
    logger.info(f"Best CV accuracy: {search.best_score_:.3f}")

    # Cross-validated predictions for ROC
    lb = LabelBinarizer()
    y_bin = lb.fit_transform(y)
    classes = lb.classes_

    y_scores_all = np.zeros_like(y_bin, dtype=float)
    y_pred_all = np.array([""] * len(y), dtype=object)

    for train_idx, test_idx in cv.split(X_scaled, y):
        X_train, X_test = X_scaled.iloc[train_idx], X_scaled.iloc[test_idx]
        y_train = y.iloc[train_idx]

        model = RandomForestClassifier(**search.best_params_, random_state=seed, n_jobs=-1)
        model.fit(X_train, y_train)

        y_scores_all[test_idx] = model.predict_proba(X_test)
        y_pred_all[test_idx] = model.predict(X_test)

    return {
        "model": best_model,
        "X": X_scaled,
        "y": y,
        "y_bin": y_bin,
        "y_scores": y_scores_all,
        "y_pred": y_pred_all,
        "classes": classes,
        "feature_importances": best_model.feature_importances_,
        "feature_names": X.columns.tolist(),
        "best_score": search.best_score_,
    }


def fig5a_roc(results: dict):
    """ROC curves (one-vs-rest)."""
    logger.info("  Fig 5a: ROC curves")

    classes = results["classes"]
    y_bin = results["y_bin"]
    y_scores = results["y_scores"]

    fig, ax = plt.subplots(figsize=(8, 7))

    for i, cls in enumerate(classes):
        fpr, tpr, _ = roc_curve(y_bin[:, i], y_scores[:, i])
        auc_val = roc_auc_score(y_bin[:, i], y_scores[:, i])
        color = LME_COLORS.get(cls, None)
        ax.plot(fpr, tpr, label=f"{cls} (AUC={auc_val:.2f})", color=color, linewidth=2)

    # Macro average
    macro_auc = roc_auc_score(y_bin, y_scores, average="macro", multi_class="ovr")

    ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title(f"ROC Curves — LME Classification (Macro AUC={macro_auc:.2f})")
    ax.legend(loc="lower right", fontsize=9)
    plt.tight_layout()

    out = FIG_DIR / "fig5a_roc.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig5b_feature_importance(results: dict, top_n: int = 20):
    """Feature importance barplot."""
    logger.info("  Fig 5b: Feature importance")

    importances = results["feature_importances"]
    features = results["feature_names"]

    idx = np.argsort(importances)[::-1][:top_n]

    fig, ax = plt.subplots(figsize=(8, max(5, top_n * 0.3)))
    ax.barh(
        range(len(idx)),
        importances[idx][::-1],
        color="#4575b4",
        edgecolor="white",
    )
    ax.set_yticks(range(len(idx)))
    ax.set_yticklabels([features[i] for i in idx[::-1]], fontsize=9)
    ax.set_xlabel("Feature Importance")
    ax.set_title(f"Top {top_n} Features — Random Forest")
    plt.tight_layout()

    out = FIG_DIR / "fig5b_feature_importance.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig5c_confusion(results: dict):
    """Confusion matrix."""
    logger.info("  Fig 5c: Confusion matrix")

    classes = results["classes"]
    cm = confusion_matrix(results["y"], results["y_pred"], labels=classes)
    cm_norm = cm.astype(float) / cm.sum(axis=1, keepdims=True)

    fig, ax = plt.subplots(figsize=(8, 7))
    sns.heatmap(
        cm_norm, annot=True, fmt=".2f", cmap="Blues", ax=ax,
        xticklabels=classes, yticklabels=classes,
        cbar_kws={"label": "Proportion"},
    )
    # Add raw counts
    for i in range(len(classes)):
        for j in range(len(classes)):
            ax.text(j + 0.5, i + 0.7, f"(n={cm[i, j]})",
                    ha="center", va="center", fontsize=7, color="gray")

    ax.set_xlabel("Predicted")
    ax.set_ylabel("True")
    ax.set_title("Confusion Matrix (normalized)")
    plt.tight_layout()

    out = FIG_DIR / "fig5c_confusion.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig5d_ihc_validation(config: dict, results: dict):
    """IHC validation scatter."""
    logger.info("  Fig 5d: IHC validation")

    ihc_path = PROJECT_DIR / config["metadata"].get("ihc_337", "")
    if not ihc_path.exists():
        logger.warning(f"  IHC file not found: {ihc_path}")
        return

    ihc = pd.read_csv(ihc_path)
    logger.info(f"  IHC data: {ihc.shape}")

    # Placeholder: scatter of IHC vs model predictions
    # Actual implementation depends on IHC column structure
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.text(0.5, 0.5, "IHC Validation\n(requires IHC data alignment)",
            ha="center", va="center", transform=ax.transAxes, fontsize=12)
    ax.set_title("IHC Validation (placeholder)")
    plt.tight_layout()

    out = FIG_DIR / "fig5d_ihc_validation.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig5e_metrics_table(results: dict):
    """Classification metrics summary."""
    logger.info("  Fig 5e: Metrics table")

    report = classification_report(
        results["y"], results["y_pred"],
        output_dict=True, zero_division=0,
    )
    report_df = pd.DataFrame(report).T

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.axis("off")
    table = ax.table(
        cellText=report_df.round(3).values,
        colLabels=report_df.columns,
        rowLabels=report_df.index,
        cellLoc="center",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.2, 1.5)
    ax.set_title("Classification Report", fontsize=12, pad=20)
    plt.tight_layout()

    out = FIG_DIR / "fig5e_metrics_table.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()

    # Also save CSV
    report_df.to_csv(FIG_DIR / "fig5e_classification_report.csv")
    logger.info(f"    Saved: {out}")


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    config = load_config()

    X, y = load_data(config)
    if X is None:
        logger.error("Cannot generate Figure 5 without TME z-scores + LME labels")
        return

    results = train_and_evaluate(X, y, config)

    fig5a_roc(results)
    fig5b_feature_importance(results)
    fig5c_confusion(results)
    fig5d_ihc_validation(config, results)
    fig5e_metrics_table(results)

    logger.info("Figure 5 complete.")


if __name__ == "__main__":
    main()
