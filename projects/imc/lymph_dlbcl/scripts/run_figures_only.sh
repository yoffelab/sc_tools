#!/bin/bash
# Run figure scripts only (data already exists from prior runs)
cd /athena/elementolab/scratch/juk4007/sc_tools/projects/imc/lymph_dlbcl

LOGDIR=outputs
mkdir -p $LOGDIR

echo "=== Starting figures $(date) ==="

for f in scripts/fig{1,2,3,4,5}_*.py scripts/supp_fig{1,2,3,4,5,6,7,8}_*.py; do
    if [ -f "$f" ]; then
        base=$(basename "$f" .py)
        echo "Running: $f"
        python3 "$f" > "$LOGDIR/${base}.log" 2>&1 && echo "  OK" || echo "  FAILED (see $LOGDIR/${base}.log)"
    fi
done

echo "=== All figures complete $(date) ==="
echo "PDF count: $(find figures/manuscript -name '*.pdf' | wc -l)"
find figures/manuscript -name '*.pdf' -printf '%p (%s bytes)\n' | sort
