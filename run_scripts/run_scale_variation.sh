#!/bin/bash
# Runs run_nucleus.sh with the fragmentation scale Q = factor * m_T, for each factor in SCALE_FACTORS.
# Usage: FRAG_TYPE=BCFY ./run_scripts/run_scale_variation.sh   (output in output/<CHANNEL>/scale_variation/<FRAG>/)

set -e
cd "$(dirname "$0")/.."
source run_scripts/config.sh

OUTBASE=${OUTBASE:-$OUTPUT_ROOT/$channel_tag/scale_variation/$FRAG_TYPE}

for factor in $SCALE_FACTORS; do
	echo "=== $FRAG_TYPE, scale factor $factor ($(date)) ==="
	SCALE_FACTOR="$factor" OUTDIR="$OUTBASE/factor_${factor}" \
		./run_scripts/run_nucleus.sh
done

echo "Finished at $(date)"
