#!/bin/bash

# Usage: Y_VALS="0.0 2.0 4.0" NPT=60 CALLS=1e6 ./local_workflows/run_fixed_qp.sh
# Runs charm_fixed_qp once per rapidity, writing output/charm/charm_fixed_qp_y<ytag>.dat.

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
