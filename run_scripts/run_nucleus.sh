#!/bin/bash
# D0 cross section for every Glauber sample b_d, pD0 and y, in output/<CHANNEL>/central_values/<FRAG>/.
# Usage: CHANNEL=0n0n FRAG_TYPE=BCFY ./run_scripts/run_nucleus.sh

set -e
cd "$(dirname "$0")/.."
source run_scripts/config.sh

require_dir()  { [[ -d "$1" ]] || { echo "Error: $1 does not exist. $2" >&2; exit 1; }; }
# "-1.5" -> "-15": the y/pT tag in file names
tag() { echo "$1" | tr -d '.'; }
# Wait until fewer than CORES jobs are running.
throttle() { while (( $(jobs -r | wc -l) >= CORES )); do sleep 0.2; done; }
# Processes to run (PROCESS = exclusive, diffractive or both)
processes() { [[ "$PROCESS" == "both" ]] && echo "exclusive diffractive" || echo "$PROCESS"; }
# One column of a program's output line
data_column() { awk -v c="$2" '$1 !~ /^#/ {print $c}' <<< "$1"; }

TARGET=AA   # always the Pb+Pb photon flux
DIPOLE_DIR=${DIPOLE_DIR:-data/$NUCLEUS/mve}
OUTDIR=${OUTDIR:-$OUTPUT_ROOT/$channel_tag/central_values/$FRAG_TYPE}

require_dir "$DIPOLE_DIR" "(expected Glauber samples glauber_mve_<b_d>)"
mkdir -p "$OUTDIR"

outfile() { echo "$OUTDIR/D0_${1}_${FRAG_TYPE}_${channel_tag}_${NUCLEUS}${flux_tag}_y$(tag "$2").dat"; }

for proc in $(processes); do
	if [[ "$proc" == "exclusive" ]]; then
		proc_calls="${CALLS_EXCL} below ${EXCL_PT_THRESHOLD} GeV, ${CALLS_EXCL_FACTORIZED} above (see EXCL_PT_THRESHOLD)"
	else
		proc_calls=$CALLS_DIFF
	fi
	for y in $Y_VALS; do
		{
			echo "# ${proc} D0 cross section, ${NUCLEUS} target, ${FRAG_TYPE} fragmentation, ${channel_tag} channel"
			echo "# generated      : $(date '+%Y-%m-%d %H:%M:%S %Z')"
			echo "# flux_model     : ${FLUX_MODEL} (TARGET=AA)"
			echo "# scale_factor   : ${SCALE_FACTOR} (fragmentation scale Q = SCALE_FACTOR * m_T)"
			echo "# calls          : ${proc_calls} (VEGAS calls per point, see CALLS_EXCL/CALLS_EXCL_FACTORIZED/CALLS_DIFF)"
			echo "# dipole samples : ${DIPOLE_DIR}/glauber_mve_<b_d>"
			echo "# pD0 sweep      : $(echo $PT_VALS)"
			echo "# fixed rapidity y : ${y}"
			echo "# ============================================================"
			echo "# b_d  pD0  dsigma_dyd2pD0"
		} > "$(outfile "$proc" "$y")"
	done
done

echo "Running D0 over $DIPOLE_DIR, pD0 in {$(echo $PT_VALS)}, y in {$Y_VALS} -> $OUTDIR"
echo "frag_type=$FRAG_TYPE channel=$CHANNEL flux_model=$FLUX_MODEL scale_factor=$SCALE_FACTOR process=$PROCESS"

# D0 prints "pD0  exclusive  diffractive"
run_one_point() {
	local dfile="$1" b_d="$2" pt="$3" y="$4"
	local result excl diff
	result=$(./build/bin/D0 "$dfile" "$pt" "$y" "$FRAG_TYPE" "$CHANNEL")
	excl=$(data_column "$result" 2)
	diff=$(data_column "$result" 3)
	if [[ -z "$excl" || -z "$diff" ]]; then
		echo "Warning: D0 $dfile $pt $y produced no data line -- skipping this point." >&2
		return
	fi
	[[ "$PROCESS" != "diffractive" ]] && echo "$b_d  $pt  $excl" >> "$(outfile exclusive "$y")"
	[[ "$PROCESS" != "exclusive" ]]   && echo "$b_d  $pt  $diff" >> "$(outfile diffractive "$y")"
	return 0
}

for dfile in "$DIPOLE_DIR"/glauber_mve_*; do
	b_d=$(basename "$dfile" | sed 's/glauber_mve_//')
	for pt in $PT_VALS; do
		for y in $Y_VALS; do
			run_one_point "$dfile" "$b_d" "$pt" "$y" &
			throttle
		done
	done
done
wait

echo "Done. Next: python3 D0.py (from plotting_scripts/) to integrate over b_d and plot."
