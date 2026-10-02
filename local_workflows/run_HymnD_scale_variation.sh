#!/bin/bash
# Usage: ./local_workflows/run_HymnD_scale_variation.sh   [SCALE_FACTORS="0.5 2.0"]
# Runs run_many_nucleus.sh with Q = SCALE_FACTOR * mt0 for the central HymnD member. Output in $OUTBASE/factor_<f>/output/.

set -e

HYMND_DIR=${HYMND_DIR:-inputs/HymnD}
HYMND_SET=${HYMND_SET:-prompt-D0-1-109}
OUTBASE=${OUTBASE:-output/HymnD_scale}
SCALE_FACTORS=${SCALE_FACTORS:-"0.5 2.0"}

member_file="$HYMND_DIR/${HYMND_SET}_0000.dat"
if [[ ! -f "$member_file" ]]; then
	echo "Error: central member file not found: $member_file" >&2
	exit 1
fi

mkdir -p "$OUTBASE"

for factor in $SCALE_FACTORS; do
	echo "=== scale factor $factor ($(date)) ==="
	PT_VALS="$(seq 0.2 0.1 2.0) $(seq 2.5 0.5 12.0)" \
	FRAG_TYPE=HymnD \
	HYMND_FILE="$member_file" \
	SCALE_FACTOR="$factor" \
	Y_VALS="-1.0 0.0 1.0 2.0 3.0" \
	PROCESS=both \
	CALLS_EXCL=1e5 CALLS_EXCL_FACTORIZED=2e3 CALLS_DIFF=1e5 \
	OUTDIR="$OUTBASE/factor_${factor}" \
	./local_workflows/run_many_nucleus.sh
done

echo "Finished at $(date)"
