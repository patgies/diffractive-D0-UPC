#!/bin/bash
# D0 cross section on a proton. TARGET=AA (default): R_pA baseline with the Pb+Pb flux,
# in output/<CHANNEL>/proton_baseline/<FRAG>/. TARGET=pA: p+Pb at 8.16 TeV, in output/pPb/<FRAG>/.

# Defaults for p+Pb
if [[ "${TARGET:-AA}" == "pA" ]]; then
	: "${FLUX_MODEL:=WS}"
	: "${SIGMA_NN:=99}"
	export SIGMA_NN
fi
set -e
cd "$(dirname "$0")/.."
source run_scripts/config.sh

require_file() { [[ -f "$1" ]] || { echo "Error: $1 does not exist. $2" >&2; exit 1; }; }
# "-1.5" -> "-15": the y/pT tag in file names
tag() { echo "$1" | tr -d '.'; }
# Wait until fewer than CORES jobs are running.
throttle() { while (( $(jobs -r | wc -l) >= CORES )); do sleep 0.2; done; }
# Processes to run (PROCESS = exclusive, diffractive or both)
processes() { [[ "$PROCESS" == "both" ]] && echo "exclusive diffractive" || echo "$PROCESS"; }
# One column of a program's output line
data_column() { awk -v c="$2" '$1 !~ /^#/ {print $c}' <<< "$1"; }

DIPOLE_FILE=${DIPOLE_FILE:-data/proton/mve.dat}
if [[ "$TARGET" == "pA" ]]; then
	OUTDIR=${OUTDIR:-$OUTPUT_ROOT/pPb/$FRAG_TYPE}
	prefix=D0_pA
	name_tag=""                     # no neutron class in p+Pb
	channel_text=""
	description="p+Pb UPC (TARGET=pA): lead emits, proton target, 8.16 TeV, Gamma_pA [sigma_NN=${SIGMA_NN} mb], flux_model=${FLUX_MODEL}, no EMD"
else
	OUTDIR=${OUTDIR:-$OUTPUT_ROOT/$channel_tag/proton_baseline/$FRAG_TYPE}
	prefix=D0_proton_baseline
	name_tag="_${channel_tag}"
	channel_text=", ${channel_tag} channel"
	description="proton-target baseline for R_pA (same nuclear photon flux as Pb+Pb, TARGET=AA, flux_model=${FLUX_MODEL})"
fi

require_file "$DIPOLE_FILE"
mkdir -p "$OUTDIR"

outfile() { echo "$OUTDIR/${prefix}_${1}_${FRAG_TYPE}${name_tag}_y$(tag "$2").dat"; }

for proc in $(processes); do
	for y in $Y_VALS; do
		{
			echo "# ${description}"
			echo "# ${proc} D0 cross section, proton target, ${FRAG_TYPE} fragmentation${channel_text}"
			echo "# generated      : $(date '+%Y-%m-%d %H:%M:%S %Z')"
			echo "# scale_factor   : ${SCALE_FACTOR} (fragmentation scale Q = SCALE_FACTOR * m_T)"
			echo "# dipole file    : ${DIPOLE_FILE}"
			echo "# pD0 sweep      : $(echo $PT_VALS)"
			echo "# fixed rapidity y : ${y}"
			echo "# ============================================================"
			echo "# pD0  dsigma_dyd2pD0_raw"
		} > "$(outfile "$proc" "$y")"
	done
done

echo "Running D0 (proton target, TARGET=$TARGET) over pD0 in {$(echo $PT_VALS)}, y in {$Y_VALS} -> $OUTDIR"

run_one_point() {
	local pt="$1" y="$2"
	local result excl diff
	result=$(./build/bin/D0 "$DIPOLE_FILE" "$pt" "$y" "$FRAG_TYPE" "$CHANNEL")
	excl=$(data_column "$result" 2)
	diff=$(data_column "$result" 3)
	if [[ -z "$excl" || -z "$diff" ]]; then
		echo "Warning: D0 $pt $y produced no data line -- skipping this point." >&2
		return
	fi
	[[ "$PROCESS" != "diffractive" ]] && echo "$pt  $excl" >> "$(outfile exclusive "$y")"
	[[ "$PROCESS" != "exclusive" ]]   && echo "$pt  $diff" >> "$(outfile diffractive "$y")"
	return 0
}

for y in $Y_VALS; do
	for pt in $PT_VALS; do
		run_one_point "$pt" "$y" &
		throttle
	done
done
wait

[[ "$TARGET" == "AA" ]] && echo "Done. Next: python3 RpA.py (from plotting_scripts/) to compute R_pA and plot." || echo "Done."
