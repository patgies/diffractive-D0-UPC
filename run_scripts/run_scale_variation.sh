#!/bin/bash
# Runs the calculation with the fragmentation scale Q = factor * m_T, for each factor in SCALE_FACTORS.
# Pb+Pb (run_nucleus.sh), in output/<CHANNEL>/scale_variation/<FRAG>/factor_<factor>/:
#   FRAG_TYPE=BCFY ./run_scripts/run_scale_variation.sh
# p+Pb (run_proton.sh), in output/pPb/scale_variation/<FRAG>/factor_<factor>/:
#   TARGET=pA FRAG_TYPE=BCFY ./run_scripts/run_scale_variation.sh

# Defaults for p+Pb (the same as in run_proton.sh)
if [[ "${TARGET:-AA}" == "pA" ]]; then
	: "${FLUX_MODEL:=WS}"
	: "${SIGMA_NN:=99}"
	export SIGMA_NN
fi
set -e
cd "$(dirname "$0")/.."
source run_scripts/config.sh

if [[ "$TARGET" == "pA" ]]; then
	OUTBASE=${OUTBASE:-$OUTPUT_ROOT/pPb/scale_variation/$FRAG_TYPE}
	script=./run_scripts/run_proton.sh
else
	OUTBASE=${OUTBASE:-$OUTPUT_ROOT/$channel_tag/scale_variation/$FRAG_TYPE}
	script=./run_scripts/run_nucleus.sh
fi

for factor in $SCALE_FACTORS; do
	echo "=== $FRAG_TYPE, scale factor $factor ($(date)) ==="
	SCALE_FACTOR="$factor" OUTDIR="$OUTBASE/factor_${factor}" "$script"
done

echo "Finished at $(date)"
