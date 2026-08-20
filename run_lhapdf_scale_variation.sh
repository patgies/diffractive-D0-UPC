#!/bin/bash
# Usage: ./run_lhapdf_scale_variation.sh
#
# Runs run_many_nucleus.sh twice with the LHAPDF fragmentation scale
# Q = SCALE_FACTOR * mt0 (mt0 = sqrt(pD0^2 + m_charm^2)) set to the
# conventional up/down variation (0.5x and 2x) around the central scale
# (SCALE_FACTOR=1.0, already produced by the normal FRAG_TYPE=LHAPDF run
# into files/), using only the central LHAPDF member (member 0000 -- this
# is a scale-convention envelope, not a PDF/FF-fit uncertainty, so unlike
# run_lhapdf_members_roihu.sbatch it does not need to be redone per replica
# member).
#
# Output lands in $OUTBASE/factor_<0.5|2.0>/files/D0_<process>_LHAPDF_<channel>_<NUCLEUS>_y<Y>.dat
#
# cross_section.py combines these with the central (SCALE_FACTOR=1.0) run
# to build a min/max scale-variation envelope, separate from the replica
# (mean +/- std) band built from run_lhapdf_members_roihu.sbatch's output.
#
#   SCALE_FACTORS="0.5 2.0" ./run_lhapdf_scale_variation.sh   # default
#   SCALE_FACTORS="0.25 4.0" ./run_lhapdf_scale_variation.sh  # wider test

set -e

LHAPDF_DIR=${LHAPDF_DIR:-data/prompt-D0-1-109}
LHAPDF_SET=${LHAPDF_SET:-prompt-D0-1-109}
OUTBASE=${OUTBASE:-files/lhapdf_scale}
SCALE_FACTORS=${SCALE_FACTORS:-"0.5 2.0"}

member_file="$LHAPDF_DIR/${LHAPDF_SET}_0000.dat"
if [[ ! -f "$member_file" ]]; then
	echo "Error: central member file not found: $member_file" >&2
	exit 1
fi

mkdir -p "$OUTBASE"

for factor in $SCALE_FACTORS; do
	echo "=== scale factor $factor ($(date)) ==="
	PT_VALS="$(seq 0.2 0.1 2.0) $(seq 2.5 0.5 12.0)" \
	FRAG_TYPE=LHAPDF \
	LHAPDF_FILE="$member_file" \
	SCALE_FACTOR="$factor" \
	Y_VALS="-1.0 0.0 1.0 2.0 3.0" \
	PROCESS=both \
	CALLS_EXCL=1e5 CALLS_EXCL_FACTORIZED=2e3 CALLS_DIFF=1e5 \
	OUTDIR="$OUTBASE/factor_${factor}" \
	./run_many_nucleus.sh
done

echo "Finished at $(date)"
