#!/bin/bash

# Usage: ./run_proton_baseline.sh
#        Y_VALS="0.0 2.0" PT_VALS="1.0 4.0" ./run_proton_baseline.sh
#
# Computes the R_pA nuclear-modification-factor baseline: dsigma_pA(B_perp >
# 2R), i.e. exclusive+diffractive D0 production off a BARE PROTON, but using
# the SAME nuclear photon flux as the Pb+Pb calculation (Z=82, bmin=2R_Pb,
# bmax, S -- all hardcoded in src/main.cpp regardless of which dipole file
# is given). This is what the R_pA formula's "not directly measurable"
# baseline means: same photon source geometry as Pb+Pb, but a proton
# target instead of the nucleus.
#
# Unlike run_many_nucleus.sh, there's no Glauber b-average here: a single
# proton has no ensemble of nucleon positions to sample over (the internal
# photon-flux b-integral, bmin to bmax, already happens inside each D0
# call regardless of target). So this loops once per (pD0, y) point,
# not once per (Glauber sample, pD0, y) -- much cheaper than the Pb+Pb sweep.
#
# Env vars:
#   Y_VALS   rapidities to scan (default matches cross_section_sum.py's Y_TO_PLOT)
#   PT_VALS  pD0 values to scan, GeV (default matches run_many_nucleus.sh)
#   FRAG_TYPE  BCFY | KniehlKramer | LHAPDF (default LHAPDF)
#   CHANNEL    An0n (default) | Xn0n | PL(AnAn)
#   CORES      parallel D0 invocations (default nproc/2)
#   CALLS_EXCL, CALLS_EXCL_FACTORIZED, CALLS_DIFF, EXCL_PT_THRESHOLD (see run_many_nucleus.sh)
#   LHAPDF_FILE

set -e

Y_VALS=${Y_VALS:-"-1.0 0.0 1.0 2.0 3.0"}
PT_VALS=${PT_VALS:-"$(seq 0.2 0.1 2.0) $(seq 2.5 0.5 12.0)"}
CORES=${CORES:-$(( $(nproc) / 2 ))}
DIPOLE_FILE="data/proton/mve.dat"

frag_tag=${FRAG_TYPE:-LHAPDF}
channel=${CHANNEL:-An0n}
channel_tag=$(echo "$channel" | tr -d '() ')

CALLS_EXCL=${CALLS_EXCL:-1e5}
CALLS_EXCL_FACTORIZED=${CALLS_EXCL_FACTORIZED:-1e3}
EXCL_PT_THRESHOLD=${EXCL_PT_THRESHOLD:-3.0}
CALLS_DIFF=${CALLS_DIFF:-1e5}

export CALLS_EXCL CALLS_EXCL_FACTORIZED EXCL_PT_THRESHOLD CALLS_DIFF LHAPDF_FILE

if [[ ! -f "$DIPOLE_FILE" ]]; then
	echo "Error: $DIPOLE_FILE does not exist." >&2
	exit 1
fi

echo "Building..."
mkdir -p build
cmake -S . -B build > /dev/null
cmake --build build -j"$(nproc)" --target D0
echo "Build OK."

mkdir -p files

for proc in exclusive diffractive; do
	for y in $Y_VALS; do
		ytag=$(echo "$y" | tr -d '.')
		outfile="files/D0_proton_baseline_${proc}_${frag_tag}_${channel_tag}_y${ytag}.dat"
		{
			echo "# proton-target baseline for R_pA (same nuclear photon flux as Pb+Pb, Z=82, bmin=2R_Pb)"
			echo "# ${proc} D0 cross section, proton target, ${frag_tag} fragmentation, ${channel_tag} channel"
			echo "# generated      : $(date '+%Y-%m-%d %H:%M:%S %Z')"
			echo "# dipole file    : ${DIPOLE_FILE}"
			echo "# pD0 sweep      : $(echo $PT_VALS)"
			echo "# fixed rapidity y : ${y}"
			echo "# ============================================================"
			echo "# pD0  dsigma_dyd2pD0_raw"
		} > "$outfile"
	done
done

echo "Running D0 (proton baseline) over pD0 in {$PT_VALS}, y in {$Y_VALS} ..."

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
	echo "$pt  $excl" >> "files/D0_proton_baseline_exclusive_${frag_tag}_${channel_tag}_y${ytag}.dat"
	echo "$pt  $diff" >> "files/D0_proton_baseline_diffractive_${frag_tag}_${channel_tag}_y${ytag}.dat"
}

for y in $Y_VALS; do
	ytag=$(echo "$y" | tr -d '.')
	for pt in $PT_VALS; do
		run_one_point "$pt" "$y" "$ytag" &
		while (( $(jobs -r | wc -l) >= CORES )); do sleep 0.2; done
	done
done
wait

echo "Done. Next: python3 python/nuclear_modification_factor.py to compute R_pA and plot."
