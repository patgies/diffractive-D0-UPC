#!/bin/bash
# Photon flux scans for plotting_scripts/flux_comparison.py, in output/flux_scan/ (sigma_NN of the Starlight tables by default).
# Usage: ./run_scripts/run_flux_scan.sh   (GAMMA_AA_FILE=input/WS_photon_flux/Gamma_AA.dat for sigma_NN = 92 mb)

# Defaults of this script
: "${GAMMA_AA_FILE:=input/WS_photon_flux/Gamma_AA_sigma90.85.dat}"
: "${FLUX_CHANNELS:=AnAn An0n Xn0n 0n0n}"
set -e
cd "$(dirname "$0")/.."
source run_scripts/config.sh

require_file() { [[ -f "$1" ]] || { echo "Error: $1 does not exist. $2" >&2; exit 1; }; }
export GAMMA_AA_FILE

OUTDIR=${OUTDIR:-$OUTPUT_ROOT/flux_scan}
require_file "$GAMMA_AA_FILE"
mkdir -p "$OUTDIR"

echo "Gamma_AA table: $GAMMA_AA_FILE -> $OUTDIR"
./build/bin/scan_flux PL "PL(AnAn)" > "$OUTDIR/flux_scan_PL.dat" &
./build/bin/scan_flux EFF $FLUX_CHANNELS > "$OUTDIR/flux_scan_EFF.dat"
wait

echo "Done. Next: python3 flux_comparison.py (from plotting_scripts/) to plot."
