#!/bin/bash

# Usage: ./local_workflows/run_fixed_qp.sh
#        Y_VALS="0.0 2.0 4.0" NPT=60 CALLS=1e6 ./local_workflows/run_fixed_qp.sh
#
# Loops charm_fixed_qp once per rapidity in Y_VALS (each call already grids over
# k_perp,c internally via OpenMP -- see src/main_fixed_qp.cpp), writing
# output/charm/charm_fixed_qp_y<ytag>.dat: fixed-q+ (no photon flux), q+=2p+, exclusive
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
cmake --build build -j"$(nproc)" --target charm_fixed_qp
echo "Build OK."

mkdir -p output/charm

for y in $Y_VALS; do
	ytag=$(echo "$y" | tr -d '.')
	outfile="output/charm/charm_fixed_qp_y${ytag}.dat"
	echo "Running y=$y -> $outfile ($NPT pts, CALLS=$CALLS) ..."
	./build/bin/charm_fixed_qp "$y" "$NPT" "$CALLS" > "$outfile"
done

echo "Done. Next: python3 python/charm_fixed_qp.py to plot."
