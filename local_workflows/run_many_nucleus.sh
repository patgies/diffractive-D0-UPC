#!/bin/bash



NUCLEUS=${NUCLEUS:-Pb}
Y_VALS=${Y_VALS:-"0.0 1.0 2.0 3.0 4.0"}
PT_MIN=${PT_MIN:-0.2}
PT_STEP=${PT_STEP:-0.5}
PT_MAX=${PT_MAX:-12.0}
# PT_VALS: list of pD0 points. If not set, use the grid PT_MIN, PT_STEP, PT_MAX.
PT_VALS=${PT_VALS:-$(seq "$PT_MIN" "$PT_STEP" "$PT_MAX")}
CORES=${CORES:-$(( $(nproc) / 2 ))}
DIPOLE_DIR=${DIPOLE_DIR:-data/$NUCLEUS/mve}
# Output files go to OUTDIR/output/ (default: the current folder).
OUTDIR=${OUTDIR:-.}

frag_tag=${FRAG_TYPE:-BCFY}
channel=${CHANNEL:-An0n}
channel_tag=$(echo "$channel" | tr -d '() ')

# PROCESS=exclusive|diffractive|both: which process to run.
PROCESS=${PROCESS:-exclusive}
if [[ "$PROCESS" == "both" ]]; then
	procs="exclusive diffractive"
else
	procs="$PROCESS"
fi

# Number of VEGAS calls: CALLS_EXCL and CALLS_DIFF use CALLS if not set. From EXCL_PT_THRESHOLD
# the exclusive part uses the 3D version (CALLS_EXCL_FACTORIZED calls).

CALLS=${CALLS:-1e5}
CALLS_EXCL=${CALLS_EXCL:-1e5}
CALLS_EXCL_FACTORIZED=${CALLS_EXCL_FACTORIZED:-1e3}
EXCL_PT_THRESHOLD=${EXCL_PT_THRESHOLD:-4.0}
CALLS_DIFF=${CALLS_DIFF:-$CALLS}

export CALLS CALLS_EXCL CALLS_EXCL_FACTORIZED EXCL_PT_THRESHOLD CALLS_DIFF PROCESS HYMND_FILE

# Photon flux: EFF (default) | PL | WS. Use the same as in run_proton_baseline.sh so the flux
# cancels in R_pA. TARGET is always AA.
FLUX_MODEL=${FLUX_MODEL:-EFF}
export TARGET=AA FLUX_MODEL
# tag in the file names: nothing for the default EFF flux,
# "_STARLIGHT" etc. for the others, so runs with different fluxes keep separate files
flux_tag=$([[ "$FLUX_MODEL" == "EFF" ]] && echo "" || echo "_${FLUX_MODEL}")


mkdir -p "$OUTDIR/output"


for proc in $procs; do
	for y in $Y_VALS; do
		ytag=$(echo "$y" | tr -d '.')
		outfile="$OUTDIR/output/D0_${proc}_${frag_tag}_${channel_tag}_${NUCLEUS}${flux_tag}_y${ytag}.dat"
		if [[ "$proc" == "exclusive" ]]; then
			proc_calls="${CALLS_EXCL} below ${EXCL_PT_THRESHOLD} GeV, ${CALLS_EXCL_FACTORIZED} above (see EXCL_PT_THRESHOLD)"
		else
			proc_calls=$CALLS_DIFF
		fi
		{
			echo "# ${proc} D0 cross section, ${NUCLEUS} target, ${frag_tag} fragmentation, ${channel_tag} channel"
			echo "# generated      : $(date '+%Y-%m-%d %H:%M:%S %Z')"
			echo "# flux_model     : ${FLUX_MODEL} (TARGET=AA)"
			echo "# calls         : ${proc_calls} (VEGAS calls per point, see CALLS_EXCL/CALLS_EXCL_FACTORIZED/CALLS_DIFF)"
			echo "# dipole samples : ${DIPOLE_DIR}/glauber_mve_<b>"
			echo "# pD0 sweep      : $(echo $PT_VALS)"
			echo "# fixed rapidity y : ${y}"
			echo "# ============================================================"
			echo "# b  pD0  dsigma_dyd2pD0"
		} > "$outfile"
	done
done

echo "Running D0 over $DIPOLE_DIR, pD0 in [$PT_MIN,$PT_MAX] step $PT_STEP, y in {$Y_VALS} ..."
echo "frag_type=$frag_tag channel=$channel flux_model=$FLUX_MODEL"

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
	# Only write the files of the process that was run.
	if [[ "$PROCESS" == "both" || "$PROCESS" == "exclusive" ]]; then
		echo "$b  $pt  $excl" >> "$OUTDIR/output/D0_exclusive_${frag_tag}_${channel_tag}_${NUCLEUS}${flux_tag}_y${ytag}.dat"
	fi
	if [[ "$PROCESS" == "both" || "$PROCESS" == "diffractive" ]]; then
		echo "$b  $pt  $diff" >> "$OUTDIR/output/D0_diffractive_${frag_tag}_${channel_tag}_${NUCLEUS}${flux_tag}_y${ytag}.dat"
	fi
}

for dfile in "$DIPOLE_DIR"/glauber_mve_*; do
	b=$(basename "$dfile" | sed 's/glauber_mve_//')
	for pt in $PT_VALS; do
		for y in $Y_VALS; do
			ytag=$(echo "$y" | tr -d '.')
			run_one_point "$dfile" "$b" "$pt" "$y" "$ytag" &
			while (( $(jobs -r | wc -l) >= CORES )); do sleep 0.2; done
		done
	done
done
wait

echo "Done. Next: python3 scripts/combine_nucleus.py to integrate over b."
