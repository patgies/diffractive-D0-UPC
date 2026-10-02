#!/bin/bash

# Usage: ./local_workflows/run_pA.sh
#        Y_VALS="0.0 2.0" PT_VALS="1.0 4.0" ./local_workflows/run_pA.sh
#
# Computes exclusive+diffractive D0 production for a REAL p-Pb UPC: the lead
# ion emits the quasi-real photon and a proton is the hadronic target, using
# the actual pA flux geometry of Sec. 6.1 of arXiv:2606.05469 -- sqrt(s_NN) =
# 8.16 TeV, Gamma_pA(b) = exp(-sigma_NN*T_A(b)) in place of Gamma_AA, and no
# EMD factor (all set via TARGET=pA in src/main.cpp/photon_flux.cpp).
#
# This is NOT the same thing as run_proton_baseline.sh: that script computes
# the R_pA "denominator" baseline, which deliberately keeps the SAME nuclear
# photon flux as Pb+Pb (Z=82, bmin=2R_Pb, sqrt(s_NN)=5.36 TeV) with a proton
# target, purely for the R_pA ratio's definition. This script instead uses
# the physically correct pA flux/kinematics -- use it for actual p-Pb
# predictions, not for the R_pA baseline.
#
# The proton dipole amplitude (data/proton/mve.dat) is used as-is: sigma0 in
# its GBW normalization already represents the target proton's transverse
# area (see python/RpA.py's proton_prefactor()),
# which is the correct convention here too, since the target is a proton in
# both cases -- python/RpA.py's proton_prefactor()
# applies unchanged to this script's output (sigma0 included, no
# GEVSQR_TO_MB, no Glauber b-integral: raw D0 output is used directly).
#
# Env vars:
#   Y_VALS     rapidities to scan (default matches D0_sum.py's Y_TO_PLOT)
#   PT_VALS    pD0 values to scan, GeV (default matches run_many_nucleus.sh)
#   FRAG_TYPE  BCFY | KniehlKramer | HymnD (default HymnD)
#   CHANNEL    An0n (default) | Xn0n | 0n0n | PL(AnAn) -- unused for TARGET=pA (no EMD), kept only for the output filename tag
#   FLUX_MODEL WS (default, realistic Gamma_pA) | PL (sharp cutoff at bmin=1.1*R_A, comparison only)
#   SIGMA_NN   total pp cross section in mb (default 99, sqrt(s_NN)=8.16 TeV value from arXiv:2606.05469)
#   CORES      parallel D0 invocations (default nproc/2)
#   CALLS_EXCL, CALLS_EXCL_FACTORIZED, CALLS_DIFF, EXCL_PT_THRESHOLD (see run_many_nucleus.sh)
#   HYMND_FILE

set -e

Y_VALS=${Y_VALS:-"-1.0 0.0 1.0 2.0 3.0"}
PT_VALS=${PT_VALS:-"$(seq 0.2 0.1 2.0) $(seq 2.5 0.5 12.0)"}
CORES=${CORES:-$(( $(nproc) / 2 ))}
DIPOLE_FILE="data/proton/mve.dat"

frag_tag=${FRAG_TYPE:-HymnD}
channel=${CHANNEL:-An0n}
channel_tag=$(echo "$channel" | tr -d '() ')

CALLS_EXCL=${CALLS_EXCL:-1e5}
CALLS_EXCL_FACTORIZED=${CALLS_EXCL_FACTORIZED:-1e3}
EXCL_PT_THRESHOLD=${EXCL_PT_THRESHOLD:-3.0}
CALLS_DIFF=${CALLS_DIFF:-1e5}
SIGMA_NN=${SIGMA_NN:-99}
FLUX_MODEL=${FLUX_MODEL:-WS}

export CALLS_EXCL CALLS_EXCL_FACTORIZED EXCL_PT_THRESHOLD CALLS_DIFF HYMND_FILE
export TARGET=pA SIGMA_NN FLUX_MODEL

if [[ ! -f "$DIPOLE_FILE" ]]; then
	echo "Error: $DIPOLE_FILE does not exist." >&2
	exit 1
fi

echo "Building..."
mkdir -p build
cmake -S . -B build > /dev/null
cmake --build build -j"$(nproc)" --target D0
echo "Build OK."

mkdir -p output

for proc in exclusive diffractive; do
	for y in $Y_VALS; do
		ytag=$(echo "$y" | tr -d '.')
		outfile="output/D0_pA_${proc}_${frag_tag}_${channel_tag}_y${ytag}.dat"
		{
			echo "# p-Pb UPC (TARGET=pA): lead emits, proton target, sqrt(s_NN)=8.16 TeV,"
			echo "# Gamma_pA(b)=exp(-sigma_NN*T_A(b)) [sigma_NN=${SIGMA_NN} mb], flux_model=${FLUX_MODEL}, no EMD"
			echo "# ${proc} D0 cross section, ${frag_tag} fragmentation"
			echo "# generated      : $(date '+%Y-%m-%d %H:%M:%S %Z')"
			echo "# dipole file    : ${DIPOLE_FILE}"
			echo "# pD0 sweep      : $(echo $PT_VALS)"
			echo "# fixed rapidity y : ${y}"
			echo "# ============================================================"
			echo "# pD0  dsigma_dyd2pD0_raw"
		} > "$outfile"
	done
done

echo "Running D0 (TARGET=pA) over pD0 in {$PT_VALS}, y in {$Y_VALS} ..."

run_one_point() {
	local pt="$1" y="$2" ytag="$3"
	local result excl diff
	result=$(./build/bin/D0 "$DIPOLE_FILE" "$pt" "$y" "$frag_tag" "$channel")
	excl=$(awk '$1 !~ /^#/ {print $2}' <<< "$result")
	diff=$(awk '$1 !~ /^#/ {print $3}' <<< "$result")
	if [[ -z "$excl" || -z "$diff" ]]; then
		echo "Warning: D0 $pt $y produced no data line -- skipping this point." >&2
		return
	fi
	echo "$pt  $excl" >> "output/D0_pA_exclusive_${frag_tag}_${channel_tag}_y${ytag}.dat"
	echo "$pt  $diff" >> "output/D0_pA_diffractive_${frag_tag}_${channel_tag}_y${ytag}.dat"
}

for y in $Y_VALS; do
	ytag=$(echo "$y" | tr -d '.')
	for pt in $PT_VALS; do
		run_one_point "$pt" "$y" "$ytag" &
		while (( $(jobs -r | wc -l) >= CORES )); do sleep 0.2; done
	done
done
wait

echo "Done. output/D0_pA_{exclusive,diffractive}_${frag_tag}_${channel_tag}_y<Y>.dat written."
