#!/bin/bash
# D0 cross section differential in x_po, for every b_d, pD0, y and x_po, in output/<CHANNEL>/xpom/<FRAG>/.
# Usage: PT_VALS="1.0 4.0" Y_VALS="0.0 2.0" ./run_scripts/run_xpom.sh

# Defaults of this script (what plotting_scripts/xpom.py plots)
: "${FRAG_TYPE:=BCFY}"
: "${PT_VALS:=2.0}"
: "${Y_VALS:=-2.0 -1.0 0.0 1.0 2.0 3.0}"
set -e
cd "$(dirname "$0")/.."
source run_scripts/config.sh

require_dir()  { [[ -d "$1" ]] || { echo "Error: $1 does not exist. $2" >&2; exit 1; }; }
# "-1.5" -> "-15": the y/pT tag in file names
tag() { echo "$1" | tr -d '.'; }
# Wait until fewer than CORES jobs are running.
throttle() { while (( $(jobs -r | wc -l) >= CORES )); do sleep 0.2; done; }
# One column of a program's output line
data_column() { awk -v c="$2" '$1 !~ /^#/ {print $c}' <<< "$1"; }

TARGET=AA
DIPOLE_DIR=${DIPOLE_DIR:-data/$NUCLEUS/mve}
OUTDIR=${OUTDIR:-$OUTPUT_ROOT/$channel_tag/xpom/$FRAG_TYPE}

require_dir "$DIPOLE_DIR" "(expected Glauber samples glauber_mve_<b_d>)"
mkdir -p "$OUTDIR"

outfile() { echo "$OUTDIR/D0_${1}_xpom_${FRAG_TYPE}_${channel_tag}_${NUCLEUS}${flux_tag}_y$(tag "$2")_pt$(tag "$3").dat"; }

# x_po grid, evenly spaced in log
xpo_values=$(awk -v a="$XPO_MIN" -v b="$XPO_MAX" -v n="$XPO_N" 'BEGIN {
	la = log(a); lb = log(b)
	for (i = 0; i < n; i++) printf "%s%.8e", (i ? " " : ""), exp(la + (lb - la) * i / (n - 1))
	print ""
}')

for y in $Y_VALS; do
	for pt in $PT_VALS; do
		for proc in exclusive diffractive; do
			{
				echo "# ${proc} D0 dsigma/(d2pD0 dy dx_po), ${NUCLEUS} target, ${FRAG_TYPE} fragmentation, ${channel_tag} channel"
				echo "# generated      : $(date '+%Y-%m-%d %H:%M:%S %Z')"
				echo "# flux_model     : ${FLUX_MODEL} (TARGET=AA)"
				echo "# scale_factor   : ${SCALE_FACTOR} (fragmentation scale Q = SCALE_FACTOR * m_T)"
				echo "# calls          : ${CALLS_EXCL_FACTORIZED} (exclusive) / ${CALLS_DIFF} (diffractive) VEGAS calls per point"
				echo "# dipole samples : ${DIPOLE_DIR}/glauber_mve_<b_d>"
				echo "# x_po grid      : ${XPO_MIN} to ${XPO_MAX}, ${XPO_N} log-spaced points"
				echo "# fixed rapidity y : ${y}"
				echo "# fixed pD0        : ${pt}"
				echo "# ============================================================"
				echo "# b_d  x_po  dsigma_dxpo"
			} > "$(outfile "$proc" "$y" "$pt")"
		done
	done
done

echo "Running D0_xpom over $DIPOLE_DIR, pD0 in {$PT_VALS}, y in {$Y_VALS}, x_po in [$XPO_MIN,$XPO_MAX] ($XPO_N pts) -> $OUTDIR"
echo "frag_type=$FRAG_TYPE channel=$CHANNEL flux_model=$FLUX_MODEL"

# D0_xpom prints "pD0  x_po  exclusive  diffractive"
run_one_point() {
	local dfile="$1" b_d="$2" pt="$3" y="$4" xpo="$5"
	local result excl diff
	result=$(./build/bin/D0_xpom "$dfile" "$pt" "$y" "$xpo" "$FRAG_TYPE" "$CHANNEL")
	excl=$(data_column "$result" 3)
	diff=$(data_column "$result" 4)
	if [[ -z "$excl" || -z "$diff" ]]; then
		echo "Warning: D0_xpom $dfile $pt $y $xpo produced no data line -- skipping this point." >&2
		return
	fi
	echo "$b_d  $xpo  $excl" >> "$(outfile exclusive "$y" "$pt")"
	echo "$b_d  $xpo  $diff" >> "$(outfile diffractive "$y" "$pt")"
}

for dfile in "$DIPOLE_DIR"/glauber_mve_*; do
	b_d=$(basename "$dfile" | sed 's/glauber_mve_//')
	for y in $Y_VALS; do
		for pt in $PT_VALS; do
			for xpo in $xpo_values; do
				run_one_point "$dfile" "$b_d" "$pt" "$y" "$xpo" &
				throttle
			done
		done
	done
done
wait

echo "Done. Next: python3 xpom.py (from plotting_scripts/) to integrate over b_d and plot."
