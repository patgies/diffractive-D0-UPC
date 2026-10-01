#!/bin/bash

# Usage: NUCLEUS=Pb ./local_workflows/run_many_xpom.sh
#        PT_VALS="1.0 4.0" Y_VALS="0.0 2.0" ./local_workflows/run_many_xpom.sh
#
# Loops D0_xpom (one point at a time: one dipole file, one pD0, one y,
# one x_po) over every Glauber-sampled dipole file in data/<NUCLEUS>/mve/,
# every pD0 in PT_VALS, every y in Y_VALS, and every x_po in a log-spaced
# grid in [XPO_MIN, XPO_MAX] -- giving both channels' dsigma/(d2pD0 dy dx_po)
# differential, unlike run_many_nucleus.sh's D0 which integrates x_po out.
# For diffractive x_po is a genuinely independent variable; for exclusive
# x_P isn't (it's a function of q+ and y), so "fixed x_P" there picks out
# one q+ via a delta function and applies the corresponding Jacobian -- see
# integrand_exclusive_xpom in src/integrand.cpp.
#
# Env vars:
#   NUCLEUS      Pb (default) | Au -- selects data/<NUCLEUS>/mve/glauber_mve_*
#   Y_VALS       rapidities to scan (default "0.0 1.0 2.0")
#   PT_VALS      pD0 values to scan, GeV (default "1.0 2.0 4.0 8.0")
#   XPO_MIN,XPO_MAX,XPO_N  log-spaced x_po grid (default 1e-6, 0.1, 25)
#   CORES        parallel D0_xpom invocations (default nproc/2)
#   FRAG_TYPE    BCFY (default) | KniehlKramer | HymnD
#   CHANNEL      An0n (default) | Xn0n | PL(AnAn)
#   CALLS, CALLS_DIFF, CALLS_EXCL_FACTORIZED (CALLS_DIFF overrides CALLS for
#                      the diffractive integrand; the exclusive one always
#                      uses the factorized analytic-r1,r2 integrand, see
#                      exclusiveCrossSection_xpom, so CALLS_EXCL is unused
#                      here -- CALLS_EXCL_FACTORIZED is what matters)
#   HYMND_FILE
#   OUTDIR       output goes to $OUTDIR/files/ (default: the repo root)
#   SKIP_BUILD   1 = use the existing build/bin/D0_xpom (cluster jobs: build
#                once with build_roihu.sh instead of every job running cmake)

set -e

NUCLEUS=${NUCLEUS:-Pb}
Y_VALS=${Y_VALS:-"-1.0 0.0 1.0 2.0 3.0 4.0"}
PT_VALS=${PT_VALS:-"1.0 2.0 4.0 8.0"}
XPO_MIN=${XPO_MIN:-1e-6}
XPO_MAX=${XPO_MAX:-0.1}
XPO_N=${XPO_N:-25}
CORES=${CORES:-$(( $(nproc) / 2 ))}
DIPOLE_DIR=${DIPOLE_DIR:-data/$NUCLEUS/mve}
OUTDIR=${OUTDIR:-.}

frag_tag=${FRAG_TYPE:-KniehlKramer}
channel=${CHANNEL:-An0n}
channel_tag=$(echo "$channel" | tr -d '() ')

# VEGAS call counts (see src/main_xpom.cpp).

CALLS=${CALLS:-1e5}
CALLS_DIFF=${CALLS_DIFF:-$CALLS}
CALLS_EXCL_FACTORIZED=${CALLS_EXCL_FACTORIZED:-1e3}

export CALLS CALLS_DIFF CALLS_EXCL_FACTORIZED HYMND_FILE

# Photon flux: EFF (default, our effective flux) or STARLIGHT (Paakkinen's
# table, AnAn/An0n only); tagged in the file names as in run_many_nucleus.sh
FLUX_MODEL=${FLUX_MODEL:-EFF}
export FLUX_MODEL
flux_tag=$([[ "$FLUX_MODEL" == "EFF" ]] && echo "" || echo "_${FLUX_MODEL}")

if [[ ! -d "$DIPOLE_DIR" ]]; then
	echo "Error: $DIPOLE_DIR does not exist (expected Glauber samples glauber_mve_<b>)." >&2
	exit 1
fi

if [[ "${SKIP_BUILD:-0}" != "1" ]]; then
	echo "Building..."
	mkdir -p build
	cmake -S . -B build > /dev/null
	cmake --build build -j"$(nproc)" --target D0_xpom
	echo "Build OK."
fi

mkdir -p "$OUTDIR/files"

# Log-spaced x_po grid, computed once (awk rather than python/numpy, which
# isn't available on the cluster compute nodes without extra modules;
# gives the same %.8e values as np.exp(np.linspace(log, log, N)))
xpo_values=$(awk -v a="$XPO_MIN" -v b="$XPO_MAX" -v n="$XPO_N" 'BEGIN {
	la = log(a); lb = log(b)
	for (i = 0; i < n; i++) printf "%s%.8e", (i ? " " : ""), exp(la + (lb - la) * i / (n - 1))
	print ""
}')

for y in $Y_VALS; do
	ytag=$(echo "$y" | tr -d '.')
	for pt in $PT_VALS; do
		pttag=$(echo "$pt" | tr -d '.')
		for process in exclusive diffractive; do
			outfile="$OUTDIR/files/D0_${process}_xpom_${frag_tag}_${channel_tag}_${NUCLEUS}${flux_tag}_y${ytag}_pt${pttag}.dat"
			{
				echo "# ${process} D0 dsigma/(d2pD0 dy dx_po), ${NUCLEUS} target, ${frag_tag} fragmentation, ${channel_tag} channel"
				echo "# generated      : $(date '+%Y-%m-%d %H:%M:%S %Z')"
				echo "# calls          : ${CALLS_EXCL_FACTORIZED} (exclusive) / ${CALLS_DIFF} (diffractive) VEGAS calls per point"
				echo "# dipole samples : ${DIPOLE_DIR}/glauber_mve_<b>"
				echo "# flux_model     : ${FLUX_MODEL}"
				echo "# x_po grid      : ${XPO_MIN} to ${XPO_MAX}, ${XPO_N} log-spaced points"
				echo "# fixed rapidity y : ${y}"
				echo "# fixed pD0        : ${pt}"
				echo "# ============================================================"
				echo "# b  x_po  dsigma_dxpo"
			} > "$outfile"
		done
	done
done

echo "Running D0_xpom over $DIPOLE_DIR, pD0 in {$PT_VALS}, y in {$Y_VALS}, x_po in [$XPO_MIN,$XPO_MAX] ($XPO_N pts) ..."
echo "frag_type=$frag_tag channel=$channel"

# D0_xpom's stdout data line is "pD0  x_po  exclusive_dxpo  diffractive_dxpo"
# (columns 3, 4) -- see src/main_xpom.cpp.
run_one_point() {
	local dfile="$1" b="$2" pt="$3" y="$4" ytag="$5" pttag="$6" xpo="$7"
	local result excl diff
	result=$(./build/bin/D0_xpom "$dfile" "$pt" "$y" "$xpo" "$frag_tag" "$channel")
	excl=$(awk '$1 !~ /^#/ {print $3}' <<< "$result")
	diff=$(awk '$1 !~ /^#/ {print $4}' <<< "$result")
	if [[ -z "$excl" || -z "$diff" ]]; then
		echo "Warning: D0_xpom $dfile $pt $y $xpo produced no data line -- skipping this point." >&2
		return
	fi
	echo "$b  $xpo  $excl" >> "$OUTDIR/files/D0_exclusive_xpom_${frag_tag}_${channel_tag}_${NUCLEUS}${flux_tag}_y${ytag}_pt${pttag}.dat"
	echo "$b  $xpo  $diff" >> "$OUTDIR/files/D0_diffractive_xpom_${frag_tag}_${channel_tag}_${NUCLEUS}${flux_tag}_y${ytag}_pt${pttag}.dat"
}

for dfile in "$DIPOLE_DIR"/glauber_mve_*; do
	b=$(basename "$dfile" | sed 's/glauber_mve_//')
	for y in $Y_VALS; do
		ytag=$(echo "$y" | tr -d '.')
		for pt in $PT_VALS; do
			pttag=$(echo "$pt" | tr -d '.')
			for xpo in $xpo_values; do
				run_one_point "$dfile" "$b" "$pt" "$y" "$ytag" "$pttag" "$xpo" &
				while (( $(jobs -r | wc -l) >= CORES )); do sleep 0.2; done
			done
		done
	done
done
wait

echo "Done. Next: python3 python/xpom_spectrum.py to integrate over b and plot (add SUM=1 for the exclusive+diffractive sum)."
