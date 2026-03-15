---
name: k-dense-ml
tier: 3
description: "[Tier 3 -- Specialized] ML model training, evaluation, and statistical analysis patterns for single-cell/spatial transcriptomics. Use when designing scVI models, tuning clustering, evaluating embeddings, computing SHAP importance, or running statistical comparisons."
allowed-tools: Read Write Edit Bash Glob Grep
---

# K-Dense ML Reference

## Core Principle

Choose the simplest model that answers the biological question. Validate assumptions before testing, report effect sizes alongside p-values, and always set random seeds.

## scvi-tools: VAE Training

### Model Selection
- **scVI**: Unsupervised integration + dimensionality reduction (default choice)
- **scANVI**: Semi-supervised when reference labels available (>30% labeled)
- **TOTALVI**: CITE-seq (RNA + protein jointly)
- **DestVI**: Spatial deconvolution with matched scRNA-seq reference

### Setup and Training Pattern
```python
import scvi

# CRITICAL: always raw counts, never log-normalized
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",           # raw count layer
    batch_key="sample_id",    # batch correction
    categorical_covariate_keys=["donor"],
    continuous_covariate_keys=["pct_counts_mt"],
)

model = scvi.model.SCVI(
    adata,
    n_latent=30,              # 10-50; 30 is robust default
    n_layers=2,               # 1-3; 2 for most datasets
    gene_likelihood="zinb",   # zinb for scRNA; nb for Visium HD
)
model.train(max_epochs=400, early_stopping=True)

# Store latent space
adata.obsm["X_scVI"] = model.get_latent_representation()
adata.layers["scvi_normalized"] = model.get_normalized_expression(library_size=1e4)
```

### Training Gotchas
- Filter genes first: `sc.pp.filter_genes(adata, min_counts=3)`
- Use HVGs: `sc.pp.highly_variable_genes(adata, n_top_genes=2000-4000)`
- GPU: `model.train(accelerator="gpu")` -- 10-50x speedup
- Save model: `model.save("model_dir/")` -- always save, retraining is expensive
- Check convergence: train_loss should plateau; validation ELBO should not diverge

### Differential Expression
```python
de = model.differential_expression(
    groupby="cell_type", group1="TypeA", group2="TypeB",
    mode="change", delta=0.25  # minimum log-fold-change threshold
)
```

## scikit-learn: Classification and Clustering

### Clustering Evaluation (Leiden Resolution Sweep)
```python
import scanpy as sc
from sklearn.metrics import silhouette_score, adjusted_rand_score

resolutions = [0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0]
results = []
for res in resolutions:
    sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}")
    labels = adata.obs[f"leiden_{res}"]
    sil = silhouette_score(adata.obsm["X_scVI"], labels, sample_size=10000)
    n_clust = labels.nunique()
    results.append({"resolution": res, "silhouette": sil, "n_clusters": n_clust})
    if "reference_labels" in adata.obs:
        results[-1]["ARI"] = adjusted_rand_score(adata.obs["reference_labels"], labels)
```

### Classification Pipeline (Cell Type Prediction)
```python
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

pipe = Pipeline([
    ("scaler", StandardScaler()),
    ("clf", RandomForestClassifier(n_estimators=200, random_state=42))
])
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
scores = cross_val_score(pipe, X, y, cv=cv, scoring="f1_macro")
# Report: f"F1 = {scores.mean():.3f} +/- {scores.std():.3f}"
```

### Key Rules
- Always `stratify=y` in train_test_split for classification
- Always `random_state=42` for reproducibility
- Scale features for SVM/KNN/PCA; not needed for tree-based models
- Use `Pipeline` to prevent data leakage in cross-validation
- Imbalanced data: use `class_weight="balanced"` or SMOTE; evaluate with F1/AUC, not accuracy

## UMAP: Embedding Visualization

### Parameter Defaults
| Purpose | n_neighbors | min_dist | n_components | metric |
|---------|-------------|----------|--------------|--------|
| Visualization | 15 | 0.1 | 2 | euclidean |
| Clustering prep | 30 | 0.0 | 10 | euclidean |
| Global structure | 100 | 0.5 | 2 | euclidean |
| Text/embeddings | 15 | 0.1 | 2 | cosine |

### Scanpy Integration
```python
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15)
sc.tl.umap(adata, min_dist=0.3, random_state=42)
```

### Gotchas
- Always set `random_state` -- UMAP is stochastic
- Standardize input when not using scVI latent (which is already normalized)
- Separate UMAP for visualization vs clustering -- different parameters
- Do NOT interpret UMAP distances as biological distances

## SHAP: Model Interpretability

### Explainer Selection
| Model type | Explainer | Speed |
|------------|-----------|-------|
| XGBoost/RF/GBM | `shap.TreeExplainer` | Fast, exact |
| Linear/Logistic | `shap.LinearExplainer` | Instant |
| Neural network | `shap.DeepExplainer` | Moderate |
| Any black-box | `shap.KernelExplainer` | Slow (last resort) |
| Unsure | `shap.Explainer` | Auto-selects |

### Workflow Pattern
```python
import shap

explainer = shap.TreeExplainer(model)
shap_values = explainer(X_test)

# Global: which features matter most?
shap.plots.beeswarm(shap_values, max_display=15)

# Local: why this prediction?
shap.plots.waterfall(shap_values[0])

# Feature relationship
shap.plots.scatter(shap_values[:, "gene_of_interest"])
```

### Gotchas
- XGBoost explains log-odds by default, not probabilities
- Background data: 100-1000 training samples for KernelExplainer
- Large datasets: compute on subset `explainer(X_test[:1000])`

## Statistical Analysis

### Test Selection Decision Tree
```
Comparing 2 groups?
  Continuous + normal  --> Independent t-test (Welch's)
  Continuous + non-normal --> Mann-Whitney U
  Paired --> Paired t-test / Wilcoxon signed-rank
  Categorical --> Chi-square / Fisher's exact

Comparing 3+ groups?
  Normal --> One-way ANOVA --> Tukey post-hoc
  Non-normal --> Kruskal-Wallis --> Dunn post-hoc

Relationship?
  Two continuous --> Pearson (normal) / Spearman (non-normal)
  Predictors --> Linear regression (continuous) / Logistic (binary)
```

### sc_tools Statistical Standards
```python
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# Always BH-FDR correction
_, padj, _, _ = multipletests(pvals, method="fdr_bh")

# Always report effect size
from pingouin import compute_effsize
d = compute_effsize(group_a, group_b, eftype="cohen")
```

### Effect Size Benchmarks
| Test | Metric | Small | Medium | Large |
|------|--------|-------|--------|-------|
| t-test | Cohen's d | 0.20 | 0.50 | 0.80 |
| ANOVA | partial eta-sq | 0.01 | 0.06 | 0.14 |
| Correlation | r | 0.10 | 0.30 | 0.50 |

### Assumption Checks (Do Before Every Test)
- **Normality**: Shapiro-Wilk (n < 5000) or visual Q-Q; n > 30 per group is often sufficient
- **Equal variance**: Levene's test; if violated, use Welch's t-test
- **Independence**: verify experimental design, not testable statistically

### Power Analysis
```python
from statsmodels.stats.power import tt_ind_solve_power
n = tt_ind_solve_power(effect_size=0.5, alpha=0.05, power=0.80, ratio=1.0)
```

## Batch Correction Evaluation

### Metrics to Report
| Metric | Measures | Good value |
|--------|----------|------------|
| iLISI (integration LISI) | Batch mixing | Higher = better mixing |
| cLISI (cell-type LISI) | Bio conservation | Higher = better preservation |
| kBET | Batch rejection rate | Low rejection = good mixing |
| ASW (batch) | Batch separation | Close to 0 = good mixing |
| ASW (cell type) | Bio separation | Close to 1 = good preservation |

### Evaluation Pattern
```python
import scib

results = scib.metrics.metrics(
    adata, adata_int,
    batch_key="batch", label_key="cell_type",
    isolated_labels_asw_=True, silhouette_=True,
    graph_conn_=True, nmi_=True, ari_=True,
)
```

## Polars: Fast DataFrames

### When to Use Over Pandas
- Metadata tables > 1M rows
- Aggregations across many groups (e.g., per-spot stats across samples)
- ETL pipelines reading/writing Parquet

### Key Patterns
```python
import polars as pl

# Lazy evaluation -- predicate/projection pushdown
result = (
    pl.scan_parquet("metadata.parquet")
    .filter(pl.col("cell_type") == "Tumor")
    .group_by("sample_id")
    .agg(pl.col("score").mean().alias("mean_score"))
    .collect()
)

# Window functions (like pandas groupby-transform)
df.with_columns(
    group_mean=pl.col("score").mean().over("cell_type")
)
```

## PyTorch Lightning: Training Deep Models

### When Relevant
- Custom VAE architectures beyond scvi-tools defaults
- Fine-tuning pretrained models (e.g., Geneformer, scGPT)
- Multi-task learning setups

### Key Patterns
```python
import lightning as L

class MyModel(L.LightningModule):
    def __init__(self):
        super().__init__()
        self.save_hyperparameters()  # always call this

    def training_step(self, batch, batch_idx):
        loss = self.compute_loss(batch)
        self.log("train_loss", loss, prog_bar=True)
        return loss

    def configure_optimizers(self):
        return torch.optim.AdamW(self.parameters(), lr=1e-3)

trainer = L.Trainer(
    max_epochs=100, accelerator="gpu",
    callbacks=[L.callbacks.EarlyStopping("val_loss", patience=10)],
    deterministic=True,  # reproducibility
)
```

### Debugging Shortcut
```python
trainer = L.Trainer(fast_dev_run=True)  # 1 batch train+val+test
```
