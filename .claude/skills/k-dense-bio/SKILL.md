---
name: k-dense-bio
tier: 2
description: "[Tier 2 — Domain] Multi-omics and systems biology API patterns for spatial transcriptomics. Covers pathway databases, deconvolution models, GRN inference, PPI networks, RNA velocity, and reference atlas queries."
allowed-tools: Read Glob Grep
---

# Multi-Omics and Systems Biology Reference

## When to use
- Pathway enrichment (KEGG, Reactome, GO, GSEA) or PPI network queries (STRING)
- Cell type deconvolution setup (Tangram, DestVI, Stereoscope, Cell2location)
- Gene regulatory network inference (GRNBoost2, GENIE3, SCENIC)
- RNA velocity / trajectory analysis (scVelo, CellRank)
- Reference atlas queries (CZ CELLxGENE Census) or pseudobulk DE (PyDESeq2)

## Quick reference

| Topic | Rule |
|-------|------|
| Deconvolution input | Raw counts only -- never log-normalized |
| scVI/DestVI setup | `setup_anndata(adata, layer="counts", batch_key=...)` |
| GRNBoost2 | Always `if __name__ == '__main__':` guard; set `seed=` for reproducibility |
| STRING species | Human=9606, Mouse=10090; threshold 400=medium, 700=high, 900=highest |
| KEGG | Org codes: `hsa`/`mmu`; max 10 entries/call; image/json: 1 only |
| Reactome | POST gene symbols one-per-line; token valid 7 days |
| Census | Always `is_primary_data == True`; pin `census_version`; check size before load |
| PyDESeq2 | Counts = samples x genes (transpose if needed); design `~batch + condition` |
| scVelo | Needs `layers["spliced"]` + `layers["unspliced"]`; start stochastic, publish dynamical |

## Pathway enrichment

### Reactome ORA
```python
import requests
resp = requests.post("https://reactome.org/AnalysisService/identifiers/",
    headers={"Content-Type": "text/plain"}, data="\n".join(gene_list))
result = resp.json()
token = result["summary"]["token"]  # valid 7 days; reuse via GET /token/{TOKEN}
enriched = [pw for pw in result["pathways"] if pw["entities"]["fdr"] < 0.05]
```

### KEGG REST API (base: https://rest.kegg.jp)
```python
import requests
requests.get("https://rest.kegg.jp/list/pathway/hsa").text      # all human pathways
requests.get("https://rest.kegg.jp/link/genes/hsa00010").text   # genes in pathway
requests.get("https://rest.kegg.jp/link/pathway/hsa:7157").text # pathways for TP53
requests.get("https://rest.kegg.jp/conv/ncbi-geneid/hsa:7157").text  # ID conversion
```

### STRING enrichment + PPI
```python
import requests, pandas as pd, io
resp = requests.post("https://string-db.org/api/tsv/enrichment",
    data={"identifiers": "\r".join(proteins), "species": 9606})
df = pd.read_csv(io.StringIO(resp.text), sep="\t")
sig = df[df["fdr"] < 0.05]  # GO, KEGG, Pfam enrichment
# Also: /api/tsv/network (required_score=700), /api/json/ppi_enrichment
```

## Deconvolution

### Tangram
```python
import tangram as tg
tg.pp_adatas(adata_sc, adata_sp, genes="markers")
tg.map_cells_to_space(adata_sc, adata_sp, mode="clusters",
    cluster_label="cell_type", density_prior="rna_count_based",
    num_epochs=500, device="cuda:0")
tg.project_cell_annotations(tg_result, adata_sp, annotation="cell_type")
```

### DestVI (scvi-tools)
```python
import scvi
# Step 1: single-cell conditional model
scvi.model.CondSCVI.setup_anndata(adata_sc, layer="counts", labels_key="cell_type")
sc_model = scvi.model.CondSCVI(adata_sc, weight_obs=True)
sc_model.train(max_epochs=300)
# Step 2: spatial deconvolution model
scvi.model.DestVI.setup_anndata(adata_sp, layer="counts")
sp_model = scvi.model.DestVI.from_rna_model(adata_sp, sc_model)
sp_model.train(max_epochs=2500)
proportions = sp_model.get_proportions()  # DataFrame: spots x cell_types
```

## GRN inference (arboreto)

```python
from arboreto.algo import grnboost2
from arboreto.utils import load_tf_names

if __name__ == "__main__":
    tf_names = load_tf_names("human_tfs.txt")
    network = grnboost2(expression_data=expr_df, tf_names=tf_names, seed=777)
    # Output: DataFrame [TF, target, importance]
    high_conf = network[network["importance"] > network["importance"].quantile(0.95)]
```
GRNBoost2 = gradient boosting (fast, >5K cells). GENIE3 = random forest (validation). Full SCENIC: arboreto -> pySCENIC cisTarget -> AUCell scoring.

## RNA velocity (scVelo)

```python
import scvelo as scv
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
# Dynamical model (publication quality; ~10-30 min for 10K cells)
scv.tl.recover_dynamics(adata, n_jobs=4)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)
scv.tl.latent_time(adata)  # shared pseudotime
scv.tl.rank_velocity_genes(adata, groupby="leiden", min_corr=0.3)
# Key outputs: obs["latent_time"], layers["velocity"], var["fit_likelihood"]
```

## CZ CELLxGENE Census

```python
import cellxgene_census
with cellxgene_census.open_soma(census_version="2023-07-25") as census:
    meta = cellxgene_census.get_obs(census, "homo_sapiens",  # check size first
        value_filter="tissue_general == 'lung' and is_primary_data == True",
        column_names=["soma_joinid"])
    adata_ref = cellxgene_census.get_anndata(census=census,  # <100K cells
        organism="Homo sapiens",
        obs_value_filter="tissue_general == 'lung' and is_primary_data == True",
        obs_column_names=["cell_type", "disease", "donor_id"])
# For >100K cells: use census[...].axis_query() with iterative processing
```

## PyDESeq2 pseudobulk DE

```python
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
# counts_df must be samples x genes; metadata indexed by sample
dds = DeseqDataSet(counts=counts_df, metadata=metadata, design="~condition", refit_cooks=True)
dds.deseq2()
ds = DeseqStats(dds, contrast=["condition", "treated", "control"])
ds.summary()
sig = ds.results_df[(ds.results_df.padj < 0.05) & (ds.results_df.log2FoldChange.abs() > 1)]
```

## Common pitfalls

| Pitfall | Fix |
|---------|-----|
| Log-normalized data to scVI/DestVI | Use raw counts in `layer="counts"` |
| GRNBoost2 hangs | Missing `if __name__ == '__main__':` guard |
| Census memory error | Filter `is_primary_data == True`; check size; `axis_query()` for >100K |
| KEGG 404 | Verify org code + entry format (`hsa:NNNNN`) |
| STRING empty network | Lower `required_score`; verify species taxon ID |
| PyDESeq2 shape error | Counts must be samples x genes, not genes x samples |
| scVelo no velocity genes | Check `layers["unspliced"]` coverage; lower `min_shared_counts` |
| Reactome no results | Use gene symbols (not Ensembl); try `/projection/` endpoint |
