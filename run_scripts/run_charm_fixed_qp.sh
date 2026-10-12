#!/bin/bash
# Charm cross section at fixed q+ (no photon flux), one file per rapidity in output/charm/.
# Usage: Y_VALS="0.0 2.0" NPT=60 ./run_scripts/run_charm_fixed_qp.sh
# With DIPOLE_DIR=data/Pb/mve: one run per Glauber sample, in output/charm/<NUCLEUS>/b<b_d>/.
# Without DIPOLE_DIR: proton, with the dipole file DIPOLE_FILE (default data/proton/mve.dat).

: "${Y_VALS:=0.0 1.0 2.0}"
: "${NPT:=40}"
set -e
cd "$(dirname "$0")/.."
source run_scripts/config.sh

# "-1.5" -> "-15": the y/pT tag in file names
tag() { echo "$1" | tr -d '.'; }

OUTDIR=${OUTDIR:-$OUTPUT_ROOT/charm}

run_one() {   # dipole file, output folder
	mkdir -p "$2"
	for y in $Y_VALS; do
		outfile="$2/charm_fixed_qp_y$(tag "$y").dat"
		echo "Running y=$y, $1 -> $outfile ($NPT pts, CALLS=$CALLS) ..."
		./build/bin/charm_fixed_qp "$y" "$NPT" "$CALLS" "$1" > "$outfile"
	done
}

if [[ -n "$DIPOLE_DIR" ]]; then
	for dfile in "$DIPOLE_DIR"/glauber_mve_*; do
		b_d=$(basename "$dfile" | sed 's/glauber_mve_//')
		run_one "$dfile" "$OUTDIR/$NUCLEUS/b$b_d"
	done
else
	run_one "${DIPOLE_FILE:-data/proton/mve.dat}" "$OUTDIR"
fi

echo "Done. Next: python3 charm_fixed_qp.py or RpA_fixed_qp.py (from plotting_scripts/) to plot."
