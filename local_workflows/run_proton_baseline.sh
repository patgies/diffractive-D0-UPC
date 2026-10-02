#!/bin/bash

# Usage: Y_VALS="0.0 2.0" PT_VALS="1.0 4.0" ./local_workflows/run_proton_baseline.sh
# R_pA baseline: proton target with the same nuclear photon flux as Pb+Pb (TARGET=AA, same FLUX_MODEL).

set -e

Y_VALS=${Y_VALS:-"-1.0 0.0 1.0 2.0 3.0"}
PT_VALS=${PT_VALS:-"$(seq 0.2 0.1 2.0) $(seq 2.5 0.5 12.0)"}
CORES=${CORES:-$(( $(nproc) / 2 ))}
DIPOLE_FILE="data/proton/mve.dat"
OUTDIR=${OUTDIR:-.}

frag_tag=${FRAG_TYPE:-HymnD}
channel=${CHANNEL:-An0n}
channel_tag=$(echo "$channel" | tr -d '() ')

CALLS_EXCL=${CALLS_EXCL:-1e5}
CALLS_EXCL_FACTORIZED=${CALLS_EXCL_FACTORIZED:-1e3}
EXCL_PT_THRESHOLD=${EXCL_PT_THRESHOLD:-3.0}
CALLS_DIFF=${CALLS_DIFF:-1e5}

FLUX_MODEL=${FLUX_MODEL:-EFF}

export CALLS_EXCL CALLS_EXCL_FACTORIZED EXCL_PT_THRESHOLD CALLS_DIFF HYMND_FILE
export TARGET=AA FLUX_MODEL

if [[ ! -f "$DIPOLE_FILE" ]]; then
	echo "Error: $DIPOLE_FILE does not exist." >&2
	exit 1
fi

if [[ "${SKIP_BUILD:-0}" != "1" ]]; then
	echo "Building..."
	mkdir -p build
	cmake -S . -B build > /dev/null
	cmake --build build -j"$(nproc)" --target D0
	echo "Build OK."
fi

mkdir -p "$OUTDIR/output"

for proc in exclusive diffractive; do
	for y in $Y_VALS; do
		ytag=$(echo "$y" | tr -d '.')
		outfile="$OUTDIR/output/D0_proton_baseline_${proc}_${frag_tag}_${channel_tag}_y${ytag}.dat"
		{
			echo "# proton-target baseline for R_pA (same nuclear photon flux as Pb+Pb, TARGET=AA, flux_model=${FLUX_MODEL})"
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
	echo "$pt  $excl" >> "$OUTDIR/output/D0_proton_baseline_exclusive_${frag_tag}_${channel_tag}_y${ytag}.dat"
	echo "$pt  $diff" >> "$OUTDIR/output/D0_proton_baseline_diffractive_${frag_tag}_${channel_tag}_y${ytag}.dat"
}

for y in $Y_VALS; do
	ytag=$(echo "$y" | tr -d '.')
	for pt in $PT_VALS; do
		run_one_point "$pt" "$y" "$ytag" &
		while (( $(jobs -r | wc -l) >= CORES )); do sleep 0.2; done
	done
done
wait

echo "Done. Next: python3 python/RpA.py to compute R_pA and plot."
