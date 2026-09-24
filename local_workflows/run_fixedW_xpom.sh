#!/bin/bash

# Usage: ./local_workflows/run_fixedW_xpom.sh
#        Y_VALS="0.0 2.0" K_VALS="1.0 4.0" ./local_workflows/run_fixedW_xpom.sh
#
# Loops D0_fixedW_xpom (one point at a time: one y, one K, one x_po) over
# every y in Y_VALS, every K in K_VALS, and every x_po in a log-spaced grid
# in [XPO_MIN, XPO_MAX] -- the fixed-W (no photon flux, proton target, no
# fragmentation) analogue of run_many_xpom.sh's D0_xpom, giving
# dsigma_fixedW/(d2K dx_po) instead of integrating x_po out (as
# run_fixedW.sh's D0_fixedW does).
#
# Env vars:
#   Y_VALS   rapidities to scan (default "0.0 1.0 2.0 3.0 4.0", matching run_fixedW.sh)
#   K_VALS   observed quark transverse momenta to scan, GeV (default "1.0 2.0 4.0 8.0")
#   XPO_MIN,XPO_MAX,XPO_N  log-spaced x_po grid (default 1e-6, 0.1, 25, matching run_many_xpom.sh)
#   CORES    parallel D0_fixedW_xpom invocations (default nproc/2)
#   CALLS    VEGAS calls per point (default 1e5)

set -e

Y_VALS=${Y_VALS:-"0.0 1.0 2.0 3.0 4.0"}
K_VALS=${K_VALS:-"1.0 2.0 4.0 8.0"}
XPO_MIN=${XPO_MIN:-1e-6}
XPO_MAX=${XPO_MAX:-0.1}
XPO_N=${XPO_N:-25}
CORES=${CORES:-$(( $(nproc) / 2 ))}
CALLS=${CALLS:-1e5}

echo "Building..."
mkdir -p build
cmake -S . -B build > /dev/null
cmake --build build -j"$(nproc)" --target D0_fixedW_xpom
echo "Build OK."

mkdir -p files

# Log-spaced x_po grid, computed once in Python (stdlib math only -- no
# numpy dependency, since a bare `module purge` compute-node environment
# on a cluster may not have it installed).
xpo_values=$(python3 -c "
import math
lo, hi, n = math.log($XPO_MIN), math.log($XPO_MAX), $XPO_N
xs = [math.exp(lo + i*(hi-lo)/(n-1)) for i in range(n)]
print(' '.join(f'{x:.8e}' for x in xs))
")

for y in $Y_VALS; do
	ytag=$(echo "$y" | tr -d '.')
	for K in $K_VALS; do
		Ktag=$(echo "$K" | tr -d '.')
		outfile="files/D0_fixedW_xpom_y${ytag}_K${Ktag}.dat"
		{
			echo "# Fixed-W (no photon flux), x_po-differential diffractive dsigma/(d2K dx_po)"
			echo "# generated   : $(date '+%Y-%m-%d %H:%M:%S %Z')"
			echo "# calls       : ${CALLS} (VEGAS calls per point)"
			echo "# x_po grid   : ${XPO_MIN} to ${XPO_MAX}, ${XPO_N} log-spaced points"
			echo "# fixed rapidity y : ${y}"
			echo "# fixed K          : ${K}"
			echo "# ============================================================"
			echo "# x_po  exclusive_fixedW  diffractive_fixedW_dxpo"
		} > "$outfile"
	done
done

echo "Running D0_fixedW_xpom over y in {$Y_VALS}, K in {$K_VALS}, x_po in [$XPO_MIN,$XPO_MAX] ($XPO_N pts) ..."

run_one_point() {
	local y="$1" ytag="$2" K="$3" Ktag="$4" xpo="$5"
	local result line
	result=$(./build/bin/D0_fixedW_xpom "$y" "$K" "$xpo" "$CALLS")
	line=$(awk '$1 !~ /^#/' <<< "$result")
	if [[ -z "$line" ]]; then
		echo "Warning: D0_fixedW_xpom $y $K $xpo produced no data line -- skipping this point." >&2
		return
	fi
	echo "$line" >> "files/D0_fixedW_xpom_y${ytag}_K${Ktag}.dat"
}

for y in $Y_VALS; do
	ytag=$(echo "$y" | tr -d '.')
	for K in $K_VALS; do
		Ktag=$(echo "$K" | tr -d '.')
		for xpo in $xpo_values; do
			run_one_point "$y" "$ytag" "$K" "$Ktag" "$xpo" &
			while (( $(jobs -r | wc -l) >= CORES )); do sleep 0.2; done
		done
	done
done
wait

echo "Done. Next: python3 python/fixedW_xpom_spectrum.py to plot."
