# DLBCL lymph_dlbcl Scripts (Runbook)

Reproducible steps for inventory, download, and AnnData build. All paths are relative to **project root** `projects/imc/lymph_dlbcl/` unless stated otherwise.

---

## Step 1: Create the inventory listing on the remote, then download it with rsync

**Remote project root:** `/home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2`

Use a two-step flow: (1) on the remote, run `find` and **write the result to a file on the remote**; (2) on your **local** machine, use **rsync** to pull that file into the project.

---

### 1a. On the remote (ssh cayuga): create the listing file

SSH in and run one of the following. Each command writes the listing to a file in your home directory so you can pull it with rsync.

**Option A — GNU `find` (Linux; preferred: relative path + size in bytes)**

```bash
REMOTE_ROOT="/home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2"
OUTPUT="$HOME/dlbcv2_csv_listing.txt"
cd "$REMOTE_ROOT"
find . -type f -name "*.csv" -printf '%P\t%s\n' | sort > "$OUTPUT"
echo "Wrote $(wc -l < "$OUTPUT") lines to $OUTPUT"
```

**Option B — Include other formats (optional: csv, h5ad, rds, h5)**

```bash
REMOTE_ROOT="/home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2"
OUTPUT="$HOME/dlbcv2_csv_listing.txt"
cd "$REMOTE_ROOT"
find . -type f \( -name "*.csv" -o -name "*.h5ad" -o -name "*.rds" -o -name "*.h5" \) -printf '%P\t%s\n' | sort > "$OUTPUT"
echo "Wrote $(wc -l < "$OUTPUT") lines to $OUTPUT"
```

**Option C — Portable find (path-only if GNU find not available)**

```bash
REMOTE_ROOT="/home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2"
OUTPUT="$HOME/dlbcv2_csv_listing.txt"
cd "$REMOTE_ROOT"
find . -type f -name "*.csv" | sed 's|^\./||' | sort > "$OUTPUT"
echo "Wrote $(wc -l < "$OUTPUT") lines to $OUTPUT"
```

The file is written to **`$HOME/dlbcv2_csv_listing.txt`** on the remote (e.g. `/home/fs01/juk4007/dlbcv2_csv_listing.txt` if your home is `/home/fs01/juk4007`).

---

### 1b. On your local machine: download the listing with rsync

From your **local** machine, in the directory where your `sc_tools` repo lives, run:

```bash
# Replace with your actual repo path if different
cd /Users/junbumkim/Documents/sc_tools

# Pull the listing from cayuga into the project metadata folder.
# If your SSH host is "cayuga" and your remote home is the default:
rsync -av cayuga:~/dlbcv2_csv_listing.txt projects/imc/lymph_dlbcl/metadata/remote_csv_listing.txt
```

**What rsync does here:**

- **`cayuga:~/dlbcv2_csv_listing.txt`** — source: file `dlbcv2_csv_listing.txt` in your home directory on the host `cayuga` (uses your SSH config).
- **`projects/imc/lymph_dlbcl/metadata/remote_csv_listing.txt`** — destination: that path in your local repo; rsync will create `metadata/` if needed and write the file there.
- **`-av`** — archive mode (`-a`) and verbose (`-v`); one file so it’s just a copy, but this is a safe default.

**If your SSH config uses a different user or host:**

- User on remote is different, e.g. `juk4007@cayuga`:
  ```bash
  rsync -av juk4007@cayuga:~/dlbcv2_csv_listing.txt projects/imc/lymph_dlbcl/metadata/remote_csv_listing.txt
  ```
- Full path on remote (if `~` is ambiguous):
  ```bash
  rsync -av cayuga:/home/fs01/juk4007/dlbcv2_csv_listing.txt projects/imc/lymph_dlbcl/metadata/remote_csv_listing.txt
  ```

After this, **`projects/imc/lymph_dlbcl/metadata/remote_csv_listing.txt`** is your local copy of the remote inventory list. Use it as input for Step 2 (build_remote_inventory.py).

---

**Expected format of `remote_csv_listing.txt`:** One line per file: either `relative_path<TAB>size_bytes` (Options A/B) or `relative_path` only (Option C). Lines starting with `#` are ignored by the inventory script. A small sample is in `metadata/remote_csv_listing_sample.txt` for testing.

---

## Step 2: Build inventory (two-step: extract, then annotate)

The inventory is produced in two steps so that **reproducible extraction** (script) is separate from **natural-language annotation** (human or AI), which benefits from understanding notebook context.

---

### Step 2a: Extract features (run script locally)

From the **repository root** (or from `projects/imc/lymph_dlbcl/`):

```bash
cd projects/imc/lymph_dlbcl
python scripts/build_remote_inventory.py
```

**Input:** `metadata/remote_csv_listing.txt` (from Step 1; may include .csv, .h5ad, .rds, .h5)  
**Output:** `metadata/remote_file_inventory.csv` with columns:
- `relative_path`, `remote_full_path`, `size_bytes`, `parent_dir`, `filename`
- `descriptor` — short heuristic label (e.g. stromal_panel_stroma_1_preprocessing)
- `data_type` — **raw** | **normalized** | **processed** | **metadata** | unknown (heuristic)
- `analysis_stage` — processing | downstream | clinical | visualization | metadata | other (heuristic)
- `detailed_annotation` — raw snippets from notebook cells where this file is read/written (input for Step 2b)
- `notebook_refs` — which notebooks reference this file
- `needs_download` — **yes** | **no** | review (heuristic; can be refined in Step 2b or Step 3)
- **`detailed_descriptor`** — **left empty by the script**; filled in Step 2b

**Optional:** To re-run the script without losing existing annotations, pass the current inventory so the script can copy over existing `detailed_descriptor` values:
```bash
python scripts/build_remote_inventory.py --merge-descriptors metadata/remote_file_inventory.csv
```

**Quick run without notebook scan** (faster; `detailed_annotation` and `notebook_refs` will be empty):
```bash
python scripts/build_remote_inventory.py --skip-notebooks
```

---

### Step 2b: Annotate detailed_descriptor (human or AI)

The **`detailed_descriptor`** column is intended to be a short, natural-language summary of what each file is and how it is used. The script does not fill it; a human or an AI (e.g. Cursor agent) should.

**How to annotate:**
1. Open `metadata/remote_file_inventory.csv`.
2. For each row, use the **`detailed_annotation`** column (code snippets from notebooks) and **`notebook_refs`** (which notebooks use this file). Optionally open the referenced notebooks to read surrounding context.
3. In **`detailed_descriptor`**, write one or two sentences that describe:
   - What the file contains (e.g. expression matrix with cells as rows, cell IDs mapping, clinical metadata, clustering labels).
   - How it is used in the pipeline (e.g. read as main expression input, joined with cell IDs, used for survival analysis).
   - Whether it is raw vs normalized/processed, and which panel (immune/stromal) if relevant.

**Example:** For a row with `relative_path` = `stroma_1_preprocessing/stroma1_quant_merged.csv` and snippets showing `pd.read_csv(..., index_col=0)` in several notebooks:
- **detailed_descriptor:** "Stromal panel (S1) merged expression matrix; cells as rows, markers as columns. Used as main input for clustering, spatial analysis, and tumor scoring. Processed/normalized counts from prior pipeline."

**Reproducibility:** The script (Step 2a) is fully reproducible. The annotations in `detailed_descriptor` are the result of a separate review step; when the listing or notebooks change, re-run Step 2a (optionally with `--merge-descriptors` to keep existing text) and then update or re-do Step 2b as needed.

**Traceability:** The script (Step 2a) writes only to `metadata/remote_file_inventory.csv`. To keep the original script output unchanged, store the annotated version in a **separate file**: **`metadata/remote_file_inventory_annotated.csv`**. Run:

```bash
python scripts/annotate_inventory_descriptors.py
```

This reads `metadata/remote_file_inventory.csv` and writes `metadata/remote_file_inventory_annotated.csv` with the `detailed_descriptor` column filled (rule-based summaries from descriptor, data_type, analysis_stage, and detailed_annotation). Use **`remote_file_inventory_annotated.csv`** for deciding what to download and for documentation.

After Step 2b, use **`metadata/remote_file_inventory_annotated.csv`** to build `metadata/files_to_download.csv` (Step 3).

---

### RDS files: usage and notebook summary

To summarize which notebooks use the downloaded RDS files (and what each RDS is), run:

```bash
python scripts/analyze_rds_notebook_usage.py
```

**Input:** `metadata/files_to_download_rds.csv` (list of RDS paths).  
**Output:**
- **`metadata/rds_notebook_usage.csv`** — One row per RDS: `relative_path`, `inferred_content` (e.g. Stromal panel (S1): Seurat object (B-cell subset)), `read_by_notebooks`, `written_by_notebooks`, `read_snippets`, `write_snippets`, `summary`. Use for understanding what each RDS contains and which notebooks depend on it.
- **`metadata/notebook_rds_summary.csv`** — One row per notebook that uses RDS: `notebook`, `rds_read`, `rds_written`, `rds_count`. Use for understanding which notebooks to run when working with the downloaded RDS files.

---

## Step 3: Identify files to download

1. Open **`metadata/remote_file_inventory_annotated.csv`** (or `metadata/remote_file_inventory.csv` if you have not run the annotation step) and use the `data_type`, `analysis_stage`, `detailed_descriptor`, `detailed_annotation`, and `needs_download` columns to identify rows that are **required** to build AnnData:
   - **Primary:** **Normalized** expression matrix with **cells as observations** (rows) for immune and stromal panels — we run the sc_tools pipeline on normalized data.
   - **Secondary:** **Raw** expression matrix with cells as observations per panel, when we can reassemble or identify it on the remote (so we have both raw and normalized; e.g. for `adata.raw` or re-normalization).
   - Include any cell metadata or sample-ID files needed to build AnnData.
2. Create `metadata/files_to_download.csv` with at least one column: `relative_path` (or `remote_full_path`). Optional: `panel` (immune/stromal), `purpose` (e.g. normalized_expression, raw_expression, cell_metadata).

Example format:

```csv
relative_path,panel,purpose
stroma_1_preprocessing/stroma1_quant_merged.csv,stromal,expression
stroma_1_preprocessing/S1_seurat_cellid.csv,stromal,cell_ids
tcell_1_preprocessing/stroma1_quant_merged.csv,immune,expression
```

---

## Step 4: Download required files (run from local)

Use a documented method so the step is reproducible.

### Option A: rsync (from laptop, with SSH access to cayuga)

**Download all RDS files** (list in `metadata/files_to_download_rds.csv`):

```bash
# From repo root (e.g. /Users/junbumkim/Documents/sc_tools)
PROJECT="projects/imc/lymph_dlbcl"
REMOTE_ROOT="cayuga:/home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2"
DEST="$PROJECT/data/downloaded"
mkdir -p "$DEST"

tail -n +2 "$PROJECT/metadata/files_to_download_rds.csv" | cut -d',' -f1 | while read -r rel; do
  [ -z "$rel" ] && continue
  mkdir -p "$DEST/$(dirname "$rel")"
  rsync -av "$REMOTE_ROOT/$rel" "$DEST/$rel"
done
```

**Or use a generic list** (`metadata/files_to_download.csv`) with column `relative_path`:

```bash
tail -n +2 "$PROJECT/metadata/files_to_download.csv" | cut -d',' -f1 | while read -r rel; do
  [ -z "$rel" ] && continue
  mkdir -p "$DEST/$(dirname "$rel")"
  rsync -av "$REMOTE_ROOT/$rel" "$DEST/$rel"
done
```

### Option B: scp (if rsync not available)

Same as above but with `scp cayuga:$REMOTE_ROOT/$rel $DEST/$rel` for each line.

### Record download manifest (traceability)

After downloading, record what was copied. From project root:

```bash
echo "remote_path,local_path,date" > metadata/download_manifest.csv
tail -n +2 metadata/files_to_download.csv | cut -d',' -f1 | while read -r rel; do
  [ -z "$rel" ] && continue
  echo "/home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2/$rel,data/downloaded/$rel,$(date -I)" >> metadata/download_manifest.csv
done
```

(Or use a small script that reads `files_to_download.csv`, checks existing files under `data/downloaded/`, and writes `download_manifest.csv` with date.)

---

## Step 5: Build AnnData (after download)

Scripts to build `results/adata.immune.*.h5ad` and `results/adata.stromal.*.h5ad` from the downloaded CSVs will be added in a follow-up. Use the inventory and notebook logic to map which CSV is expression vs metadata vs cell IDs.

---

## Traceability

- **Listing:** Keep `metadata/remote_csv_listing.txt` as the raw input (or document where it came from and the exact remote command in Journal.md).
- **Inventory:** Two-step. (1) `scripts/build_remote_inventory.py` produces `metadata/remote_file_inventory.csv` with extracted features (reproducible). (2) Annotated version with `detailed_descriptor` filled is stored in **`metadata/remote_file_inventory_annotated.csv`** (via `scripts/annotate_inventory_descriptors.py` or human/AI edit). Original script output is unchanged for traceability.
- **Download:** `metadata/files_to_download.csv` (human-curated) and `metadata/download_manifest.csv` (after download) record what was chosen and what was copied.

Log material decisions (which files were marked required, download date) in `Journal.md`.
