#!/bin/bash



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

# PROCESS=exclusive|diffractive reruns just that one process (e.g. with a
# bumped CALLS_EXCL) without recomputing -- or touching the output file of --
# the other, already-converged one. Default "both" is the original behavior.
PROCESS=${PROCESS:-exclusive}
if [[ "$PROCESS" == "both" ]]; then
	procs="exclusive diffractive"
else
	procs="$PROCESS"
fi

# VEGAS call counts (see src/main.cpp): CALLS_EXCL/CALLS_DIFF each fall back
# to CALLS if unset, so a bare CALLS=1e6 still applies to both as before.
CALLS=${CALLS:-1e5}
CALLS_EXCL=${CALLS_EXCL:-5e7}
CALLS_DIFF=${CALLS_DIFF:-$CALLS}

export CALLS CALLS_EXCL CALLS_DIFF PROCESS LHAPDF_FILE


mkdir -p files


for proc in $procs; do
	for y in $Y_VALS; do
		ytag=$(echo "$y" | tr -d '.')
		outfile="files/D0_${proc}_${frag_tag}_${channel_tag}_${NUCLEUS}_y${ytag}.dat"
		if [[ "$proc" == "exclusive" ]]; then proc_calls=$CALLS_EXCL; else proc_calls=$CALLS_DIFF; fi
		{
			echo "# ${proc} D0 cross section, ${NUCLEUS} target, ${frag_tag} fragmentation, ${channel_tag} channel"
			echo "# generated      : $(date '+%Y-%m-%d %H:%M:%S %Z')"
			echo "# calls          : ${proc_calls} (VEGAS calls per point, see CALLS_EXCL/CALLS_DIFF)"
			echo "# dipole samples : ${DIPOLE_DIR}/glauber_mve_<b>"
			echo "# pD0 sweep      : ${PT_MIN} to ${PT_MAX} GeV, step ${PT_STEP}"
			echo "# fixed rapidity y : ${y}"
			echo "# ============================================================"
			echo "# b  pD0  dsigma_dyd2pD0"
		} > "$outfile"
	done
done

echo "Running D0 over $DIPOLE_DIR, pD0 in [$PT_MIN,$PT_MAX] step $PT_STEP, y in {$Y_VALS} ..."
echo "frag_type=$frag_tag channel=$channel"

# Column 2 of D0's data line is "exclusive", column 3 is "diffractive".
run_one_point() {
	local dfile="$1" b="$2" pt="$3" y="$4" ytag="$5"
	local result excl diff
	result=$(./build/bin/D0 "$dfile" "$pt" "$y" "$frag_tag" "$channel")
	excl=$(awk '$1 !~ /^#/ {print $2}' <<< "$result")
	diff=$(awk '$1 !~ /^#/ {print $3}' <<< "$result")
	if [[ -z "$excl" || -z "$diff" ]]; then
		echo "Warning: D0 $dfile $pt $y produced no data line -- skipping this point." >&2
		return
	fi
	# Only write the file(s) for the process(es) actually computed this run --
	# PROCESS=exclusive gives diff=0 (see main.cpp), which must not clobber an
	# existing good diffractive file.
	if [[ "$PROCESS" == "both" || "$PROCESS" == "exclusive" ]]; then
		echo "$b  $pt  $excl" >> "files/D0_exclusive_${frag_tag}_${channel_tag}_${NUCLEUS}_y${ytag}.dat"
	fi
	if [[ "$PROCESS" == "both" || "$PROCESS" == "diffractive" ]]; then
		echo "$b  $pt  $diff" >> "files/D0_diffractive_${frag_tag}_${channel_tag}_${NUCLEUS}_y${ytag}.dat"
	fi
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
