#!/bin/bash

# Usage: ./run_local.sh

# Output: files/D0_exclusive.dat and
# files/D0_diffractive.dat (columns: pD0  y  dsigma).

set -e

dipole_file="${1:-data/proton/mve.dat}"
pt_min="${2:-0.2}"
pt_step="${3:-0.5}"
pt_max="${4:-12.0}"
if [ $# -ge 5 ]; then
    y_values="${*:5}"
else
    y_values="0.0 1.0 2.0 3.0 4.0"
fi

# VEGAS call counts (see src/main.cpp): CALLS_EXCL/CALLS_DIFF each fall back
# to CALLS if unset, so a bare CALLS=1e6 still applies to both as before.
CALLS=${CALLS:-1e5}
CALLS_EXCL=${CALLS_EXCL:-$CALLS}
CALLS_DIFF=${CALLS_DIFF:-$CALLS}
# PROCESS=exclusive|diffractive reruns just that one process (e.g. with a
# bumped CALLS_EXCL) without recomputing -- or touching the output file of --
# the other, already-converged one. Default "both" is the original behavior.
PROCESS=${PROCESS:-both}
export CALLS CALLS_EXCL CALLS_DIFF PROCESS LHAPDF_FILE

mkdir -p files
excl_file="files/D0_exclusive.dat"
diff_file="files/D0_diffractive.dat"

timestamp=$(date '+%Y-%m-%d %H:%M:%S %Z')
if [[ "$PROCESS" == "both" || "$PROCESS" == "exclusive" ]]; then
{
    echo "# D0 exclusive cross section, dipole file: $dipole_file"
    echo "# generated: $timestamp"
    echo "# calls: $CALLS_EXCL (VEGAS calls per point)"
    echo "# pD0  y  dsigma"
} > "$excl_file"
fi
if [[ "$PROCESS" == "both" || "$PROCESS" == "diffractive" ]]; then
{
    echo "# D0 diffractive cross section, dipole file: $dipole_file"
    echo "# generated: $timestamp"
    echo "# calls: $CALLS_DIFF (VEGAS calls per point)"
    echo "# pD0  y  dsigma"
} > "$diff_file"
fi

echo "Running dipole_file=${dipole_file}, pD0 in [${pt_min},${pt_max}] step ${pt_step}, y in {${y_values}} ..."

for pt in $(seq "$pt_min" "$pt_step" "$pt_max"); do
    for y in $y_values; do
        result=$(./build/bin/D0 "$dipole_file" "$pt" "$y")
        excl=$(awk '$1 !~ /^#/ {print $2}' <<< "$result")
        diff=$(awk '$1 !~ /^#/ {print $3}' <<< "$result")
        if [[ "$PROCESS" == "both" || "$PROCESS" == "exclusive" ]]; then
            echo "$pt  $y  $excl" >> "$excl_file"
        fi
        if [[ "$PROCESS" == "both" || "$PROCESS" == "diffractive" ]]; then
            echo "$pt  $y  $diff" >> "$diff_file"
        fi
        echo "  pD0=${pt} y=${y} -> exclusive=${excl} diffractive=${diff}"
    done
done

echo "Done"
