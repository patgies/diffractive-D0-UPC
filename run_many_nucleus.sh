#!/bin/bash

# Usage: NUCLEUS=Pb ./run_many_nucleus.sh
#        NUCLEUS=Au FRAG_TYPE=KniehlKramer ./run_many_nucleus.sh
#
# Loops d0_point (one point at a time: one dipole file, one pD0, one y)
# over every Glauber-sampled dipole file in data/<NUCLEUS>/mve/, every pD0
# in [PT_MIN, PT_MAX] (step PT_STEP), and every y in Y_VALS 
#
# Env vars:
#   NUCLEUS      Pb (default) | Au -- selects data/<NUCLEUS>/mve/glauber_mve_*
#   Y_VALS       rapidities to scan (default "0.0 1.0 2.0 3.0 4.0")
#   PT_MIN,PT_STEP,PT_MAX  pD0 sweep, GeV (default 0.2, 0.5, 12.0)
#   CORES        parallel d0_point invocations (default nproc/2)
#   FRAG_TYPE    BCFY (default) | KniehlKramer | LHAPDF
#   CHANNEL      An0n (default) | Xn0n | PL(AnAn) 
#   CALLS, LHAPDF_FILE

set -e

NUCLEUS=${NUCLEUS:-Pb}
Y_VALS=${Y_VALS:-"0.0 1.0 2.0 3.0 4.0"}
PT_MIN=${PT_MIN:-0.2}
PT_STEP=${PT_STEP:-0.5}
PT_MAX=${PT_MAX:-12.0}
CORES=${CORES:-$(( $(nproc) / 2 ))}
DIPOLE_DIR=${DIPOLE_DIR:-data/$NUCLEUS/mve}

frag_tag=${FRAG_TYPE:-BCFY}
channel=${CHANNEL:-An0n}
channel_tag=$(echo "$channel" | tr -d '() ')


export CALLS LHAPDF_FILE

if [[ ! -d "$DIPOLE_DIR" ]]; then
	echo "Error: $DIPOLE_DIR does not exist (expected Glauber samples glauber_mve_<b>)." >&2
	exit 1
fi

echo "Building..."
mkdir -p build
cmake -S . -B build > /dev/null
cmake --build build -j"$(nproc)" --target d0_point
echo "Build OK."

mkdir -p files


for proc in exclusive diffractive; do
	for y in $Y_VALS; do
		ytag=$(echo "$y" | tr -d '.')
		outfile="files/d0_point_${proc}_${frag_tag}_${channel_tag}_${NUCLEUS}_y${ytag}.dat"
		{
			echo "# ${proc} D0 cross section, ${NUCLEUS} target, ${frag_tag} fragmentation, ${channel_tag} channel"
			echo "# dipole samples : ${DIPOLE_DIR}/glauber_mve_<b>"
			echo "# pD0 sweep      : ${PT_MIN} to ${PT_MAX} GeV, step ${PT_STEP}"
			echo "# fixed rapidity y : ${y}"
			echo "# ============================================================"
			echo "# b  pD0  dsigma_dyd2pD0"
		} > "$outfile"
	done
done

echo "Running d0_point over $DIPOLE_DIR, pD0 in [$PT_MIN,$PT_MAX] step $PT_STEP, y in {$Y_VALS} ..."
echo "frag_type=$frag_tag channel=$channel"

# Column 2 of d0_point's data line is "exclusive", column 3 is "diffractive".
run_one_point() {
	local dfile="$1" b="$2" pt="$3" y="$4" ytag="$5"
	local result excl diff
	result=$(./build/bin/d0_point "$dfile" "$pt" "$y" "$frag_tag" "$channel")
	excl=$(awk '$1 !~ /^#/ {print $2}' <<< "$result")
	diff=$(awk '$1 !~ /^#/ {print $3}' <<< "$result")
	if [[ -z "$excl" || -z "$diff" ]]; then
		echo "Warning: d0_point $dfile $pt $y produced no data line -- skipping this point." >&2
		return
	fi
	echo "$b  $pt  $excl" >> "files/d0_point_exclusive_${frag_tag}_${channel_tag}_${NUCLEUS}_y${ytag}.dat"
	echo "$b  $pt  $diff" >> "files/d0_point_diffractive_${frag_tag}_${channel_tag}_${NUCLEUS}_y${ytag}.dat"
}

for dfile in "$DIPOLE_DIR"/glauber_mve_*; do
	b=$(basename "$dfile" | sed 's/glauber_mve_//')
	for pt in $(seq "$PT_MIN" "$PT_STEP" "$PT_MAX"); do
		for y in $Y_VALS; do
			ytag=$(echo "$y" | tr -d '.')
			run_one_point "$dfile" "$b" "$pt" "$y" "$ytag" &
			while (( $(jobs -r | wc -l) >= CORES )); do sleep 0.2; done
		done
	done
done
wait

echo "Done. Next: python3 scripts/combine_nucleus.py to integrate over b."
