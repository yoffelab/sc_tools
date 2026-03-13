---
type: signature
tier: knowledge
tags:
  - knowledge/signatures
---

# MSigDB Hallmark Gene Sets

50 well-defined biological state and process gene sets from MSigDB. Bundled offline in `sc_tools/data/hallmark_human.json`.

## Usage

```python
from sc_tools.tl import load_hallmark, score_signature, merge_gene_signatures

hallmark = load_hallmark()  # 50 sets, no network needed
project_sigs = load_json("metadata/gene_signatures.json")
combined = merge_gene_signatures(project_sigs, hallmark)
score_signature(adata, combined)
```

## Convention

Always apply Hallmark + project signatures together via `merge_gene_signatures()` before calling `score_signature()`. This ensures consistent column naming in `obsm['signature_score']`.

## Sets

50 sets covering: apoptosis, cell cycle, DNA repair, EMT, glycolysis, hedgehog, hypoxia, immune response, inflammation, metabolism, mTORC1, Myc, Notch, oxidative phosphorylation, p53, PI3K/Akt, reactive oxygen species, TGF-beta, TNF-alpha/NF-kB, UV response, Wnt, and more.

See `sc_tools.tl.list_gene_sets(load_hallmark())` for the full list.
