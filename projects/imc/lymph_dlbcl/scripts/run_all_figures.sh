#!/bin/bash
# Run full DLBCL figure pipeline on HPC
set -e
cd /athena/elementolab/scratch/juk4007/sc_tools/projects/imc/lymph_dlbcl

LOGDIR=outputs
mkdir -p $LOGDIR figures/manuscript

echo "=== Starting pipeline $(date) ==="

# Step 1: Build panel adata
echo "Step 1: build_panel_adata.py"
python3 scripts/build_panel_adata.py > $LOGDIR/step1_build_panel.log 2>&1
echo "  Done: $(tail -1 $LOGDIR/step1_build_panel.log)"

# Step 2: Attach clinical metadata
echo "Step 2: attach_clinical_metadata.py"
python3 scripts/attach_clinical_metadata.py > $LOGDIR/step2_attach_clinical.log 2>&1
echo "  Done: $(tail -1 $LOGDIR/step2_attach_clinical.log)"

# Step 3: Build LME classes
echo "Step 3: build_lme_classes.py"
python3 scripts/build_lme_classes.py > $LOGDIR/step3_build_lme.log 2>&1
echo "  Done: $(tail -1 $LOGDIR/step3_build_lme.log)"

echo "=== Data pipeline complete $(date) ==="

# Step 4: Run all figure scripts
for f in scripts/fig{1,2,3,4,5}_*.py scripts/supp_fig{1,2,3,4,5,6,7,8}_*.py; do
    if [ -f "$f" ]; then
        base=$(basename "$f" .py)
        echo "Running: $f"
        python3 "$f" > "$LOGDIR/${base}.log" 2>&1 && echo "  OK" || echo "  FAILED (see $LOGDIR/${base}.log)"
    fi
done

echo "=== All figures complete $(date) ==="
echo "=== PDF count: $(find figures/manuscript -name '*.pdf' | wc -l) ==="
