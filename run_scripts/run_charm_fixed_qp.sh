#!/bin/bash
# Charm cross section at fixed q+ (no photon flux), one file per rapidity in output/charm/.
# Usage: Y_VALS="0.0 2.0" NPT=60 ./run_scripts/run_charm_fixed_qp.sh

# Defaults of this script
: "${Y_VALS:=0.0 1.0 2.0}"
: "${NPT:=40}"
set -e
cd "$(dirname "$0")/.."
source run_scripts/config.sh

# "-1.5" -> "-15": the y/pT tag in file names
tag() { echo "$1" | tr -d '.'; }

OUTDIR=${OUTDIR:-$OUTPUT_ROOT/charm}
mkdir -p "$OUTDIR"

for y in $Y_VALS; do
	outfile="$OUTDIR/charm_fixed_qp_y$(tag "$y").dat"
	echo "Running y=$y -> $outfile ($NPT pts, CALLS=$CALLS) ..."
	./build/bin/charm_fixed_qp "$y" "$NPT" "$CALLS" > "$outfile"
done

echo "Done. Next: python3 charm_fixed_qp.py (from plotting_scripts/) to plot."
