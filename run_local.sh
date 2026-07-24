#!/bin/bash

# Usage: ./run_local.sh [dipole_file] [pt_min] [pt_step] [pt_max] [y_values...]

# Output: files/d0_point_exclusive.dat and
# files/d0_point_diffractive.dat (columns: pD0  y  dsigma).

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

echo "Building..."
mkdir -p build
cmake -S . -B build > /dev/null
cmake --build build -j"$(nproc)" --target d0_point
echo "Build OK."

mkdir -p files
excl_file="files/d0_point_exclusive.dat"
diff_file="files/d0_point_diffractive.dat"

{
    echo "# D0 exclusive cross section, dipole file: $dipole_file"
    echo "# pD0  y  dsigma"
} > "$excl_file"
{
    echo "# D0 diffractive cross section, dipole file: $dipole_file"
    echo "# pD0  y  dsigma"
} > "$diff_file"

echo "Running dipole_file=${dipole_file}, pD0 in [${pt_min},${pt_max}] step ${pt_step}, y in {${y_values}} ..."

for pt in $(seq "$pt_min" "$pt_step" "$pt_max"); do
    for y in $y_values; do
        result=$(./build/bin/d0_point "$dipole_file" "$pt" "$y")
        excl=$(awk '$1 !~ /^#/ {print $2}' <<< "$result")
        diff=$(awk '$1 !~ /^#/ {print $3}' <<< "$result")
        echo "$pt  $y  $excl" >> "$excl_file"
        echo "$pt  $y  $diff" >> "$diff_file"
        echo "  pD0=${pt} y=${y} -> exclusive=${excl} diffractive=${diff}"
    done
done

echo "Done. See $excl_file and $diff_file"
