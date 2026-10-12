#!/bin/bash
# Diffractive charm x_po distribution at fixed q+ (no photon flux), one file per rapidity in output/charm/xpom/.
# Usage: Y_VALS="0.0 2.0" P=2.0 NXPO=60 ./run_scripts/run_charm_fixed_qp_xpom.sh

# Defaults of this script (the rapidities of run_charm_fixed_qp.sh's plot)
: "${Y_VALS:=0.0 0.5 1.0 1.5 2.0}"
: "${P:=2.0}"
: "${NXPO:=60}"
set -e
cd "$(dirname "$0")/.."
source run_scripts/config.sh

# "-1.5" -> "-15": the y/pT tag in file names
tag() { echo "$1" | tr -d '.'; }

OUTDIR=${OUTDIR:-$OUTPUT_ROOT/charm/xpom}
mkdir -p "$OUTDIR"

for y in $Y_VALS; do
	outfile="$OUTDIR/charm_fixed_qp_xpom_y$(tag "$y")_p$(tag "$P").dat"
	echo "Running y=$y, p=$P -> $outfile ($NXPO pts, CALLS=$CALLS) ..."
	./build/bin/charm_fixed_qp_xpom "$y" "$P" "$NXPO" "$CALLS" > "$outfile"
done

echo "Done. Next: python3 charm_fixed_qp_xpom.py (from plotting_scripts/) to plot."
