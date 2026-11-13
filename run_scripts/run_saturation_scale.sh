#!/bin/bash
# Saturation scale vs x of the proton dipole and of each Pb Glauber sample, in output/saturation_scale/.
# Usage: ./run_scripts/run_saturation_scale.sh   (NS=0.5 for N(r_s) = 1/2)

: "${NS:=0.393469}"   # 1 - exp(-1/2)
set -e
cd "$(dirname "$0")/.."
source run_scripts/config.sh

OUTDIR=${OUTDIR:-$OUTPUT_ROOT/saturation_scale}
mkdir -p "$OUTDIR"

./build/bin/saturation_scale "${DIPOLE_FILE:-data/proton/mve.dat}" "$NS" > "$OUTDIR/qs_proton.dat"
for dfile in "${DIPOLE_DIR:-data/$NUCLEUS/mve}"/glauber_mve_*; do
	b_d=$(basename "$dfile" | sed 's/glauber_mve_//')
	./build/bin/saturation_scale "$dfile" "$NS" > "$OUTDIR/qs_${NUCLEUS}_b$b_d.dat"
done

echo "Done. Next: python3 correlation_limit.py (from plotting_scripts/) to plot."
