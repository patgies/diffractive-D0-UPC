#!/bin/bash

# Usage: Y_VALS="0.0 2.0" PT_VALS="1.0 4.0" ./local_workflows/run_pA.sh
# Real p-Pb UPC (TARGET=pA: 8.16 TeV, Gamma_pA, no EMD). This is not the R_pA baseline (see run_proton_baseline.sh).

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
