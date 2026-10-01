#!/bin/bash



NUCLEUS=${NUCLEUS:-Pb}
Y_VALS=${Y_VALS:-"0.0 1.0 2.0 3.0 4.0"}
PT_MIN=${PT_MIN:-0.2}
PT_STEP=${PT_STEP:-0.5}
PT_MAX=${PT_MAX:-12.0}
# PT_VALS lets you pass an explicit list of pD0 points (e.g. a finer grid at
# low pT and a coarser one at high pT) instead of one uniform PT_MIN/STEP/MAX
# sweep. If unset, falls back to the uniform seq as before.
PT_VALS=${PT_VALS:-$(seq "$PT_MIN" "$PT_STEP" "$PT_MAX")}
CORES=${CORES:-$(( $(nproc) / 2 ))}
DIPOLE_DIR=${DIPOLE_DIR:-data/$NUCLEUS/mve}
# Where output files land: OUTDIR/files/D0_<process>_..._y<Y>.dat. Defaults to
# the current directory (same as before OUTDIR existed). Useful for e.g. a
# SLURM array job giving each task its own OUTDIR so tasks don't clobber
# each other's output.
OUTDIR=${OUTDIR:-.}

frag_tag=${FRAG_TYPE:-BCFY}
channel=${CHANNEL:-An0n}
channel_tag=$(echo "$channel" | tr -d '() ')

# PROCESS=exclusive|diffractive reruns just that one process (e.g. with a
# bumped CALLS_EXCL) without recomputing  or touching the output file of 
# the other. Default "both" is the original behavior.
PROCESS=${PROCESS:-exclusive}
if [[ "$PROCESS" == "both" ]]; then
	procs="exclusive diffractive"
else
	procs="$PROCESS"
fi

# VEGAS call counts: CALLS_EXCL/CALLS_DIFF each fall back
# to CALLS if unset, so a bare CALLS=1e6 still applies to both as before.
# CALLS_EXCL used to need to be huge (5e7) because the exclusive integrand
# handed the oscillating r1, r2 Bessel integrals to VEGAS at every pt, which
# barely converged at high pt no matter how many calls were set.
# exclusiveCrossSection now switches, at EXCL_PT_THRESHOLD (default 4 GeV),
# from the 5D integrand to a factorized 3D one that solves
# r1, r2 analytically per point instead.

CALLS=${CALLS:-1e5}
CALLS_EXCL=${CALLS_EXCL:-1e5}
CALLS_EXCL_FACTORIZED=${CALLS_EXCL_FACTORIZED:-1e3}
EXCL_PT_THRESHOLD=${EXCL_PT_THRESHOLD:-4.0}
CALLS_DIFF=${CALLS_DIFF:-$CALLS}

export CALLS CALLS_EXCL CALLS_EXCL_FACTORIZED EXCL_PT_THRESHOLD CALLS_DIFF PROCESS HYMND_FILE

# Photon flux: EFF (default) | PL | WS, see src/main.cpp. Must match the
# FLUX_MODEL used for run_proton_baseline.sh, otherwise the flux no longer
# cancels in R_pA (python/nuclear_modification_factor.py). TARGET is pinned
# to AA so a stray TARGET=pA in the environment can't switch the geometry.
FLUX_MODEL=${FLUX_MODEL:-EFF}
export TARGET=AA FLUX_MODEL
# file-name tag: none for the default EFF flux (keeps the existing names),
# "_STARLIGHT" etc. otherwise, so runs with different fluxes don't overwrite
flux_tag=$([[ "$FLUX_MODEL" == "EFF" ]] && echo "" || echo "_${FLUX_MODEL}")


mkdir -p "$OUTDIR/files"


for proc in $procs; do
	for y in $Y_VALS; do
		ytag=$(echo "$y" | tr -d '.')
		outfile="$OUTDIR/files/D0_${proc}_${frag_tag}_${channel_tag}_${NUCLEUS}${flux_tag}_y${ytag}.dat"
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
	# Only write the file(s) for the process actually computed this run.
	if [[ "$PROCESS" == "both" || "$PROCESS" == "exclusive" ]]; then
		echo "$b  $pt  $excl" >> "$OUTDIR/files/D0_exclusive_${frag_tag}_${channel_tag}_${NUCLEUS}${flux_tag}_y${ytag}.dat"
	fi
	if [[ "$PROCESS" == "both" || "$PROCESS" == "diffractive" ]]; then
		echo "$b  $pt  $diff" >> "$OUTDIR/files/D0_diffractive_${frag_tag}_${channel_tag}_${NUCLEUS}${flux_tag}_y${ytag}.dat"
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
