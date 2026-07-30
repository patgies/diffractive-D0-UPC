#!/bin/bash

# Usage: ./run_fixedW.sh
#        Y_VALS="0.0 2.0 4.0" NPT=60 CALLS=1e6 ./run_fixedW.sh
#
# Loops D0_fixedW once per rapidity in Y_VALS (each call already grids over
# k_perp,c internally via OpenMP -- see src/main_fixedW.cpp), writing
# files/D0_fixedW_y<ytag>.dat: fixed-W (no photon flux), q+=3p+, exclusive
# ("missing prefactor" convention) + diffractive ("full physical cross
# section") spectra vs k_perp,c.
#
# Env vars:
#   Y_VALS   rapidities to scan (default "0.0 1.0 2.0 3.0 4.0")
#   NPT      number of k_perp,c grid points, log-spaced in [0.2,20] GeV (default 40)
#   CALLS    VEGAS calls for the diffractive integral; exclusive is deterministic (default 1e5)

set -e

Y_VALS=${Y_VALS:-"0.0 1.0 2.0 3.0 4.0"}
NPT=${NPT:-40}
CALLS=${CALLS:-1e5}

echo "Building..."
mkdir -p build
cmake -S . -B build > /dev/null
cmake --build build -j"$(nproc)" --target D0_fixedW
echo "Build OK."

mkdir -p files

for y in $Y_VALS; do
	ytag=$(echo "$y" | tr -d '.')
	outfile="files/D0_fixedW_y${ytag}.dat"
	echo "Running y=$y -> $outfile ($NPT pts, CALLS=$CALLS) ..."
	./build/bin/D0_fixedW "$y" "$NPT" "$CALLS" > "$outfile"
done

echo "Done. Next: python3 python/fixedW_spectrum.py to plot."
