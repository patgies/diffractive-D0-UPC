#!/bin/bash -l
#SBATCH --job-name=D0_lhapdf_band
#SBATCH --account=lappi
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=1000M
#SBATCH --array=0-100
#SBATCH --output=logs/lhapdf_%A_%a.out
#SBATCH --error=logs/lhapdf_%A_%a.err

# Usage (run build_roihu.sh first -- this script does NOT build):
#   sbatch run_lhapdf_members_roihu.sbatch
#   sbatch --array=0-9 run_lhapdf_members_roihu.sbatch          # quick subset test
#   sbatch --export=ALL,LHAPDF_DIR=/other/path run_lhapdf_members_roihu.sbatch
#
# One array task per LHAPDF replica member (0-100: member 0 is the central
# value, 1-100 are the replicas -- see data/prompt-D0-1-109/*.info,
# ErrorType=replicas). Each task runs the exact same diffractive+exclusive
# sweep as the local member-0000 run (fine pT grid below 2 GeV, coarser
# above, EXCL_PT_THRESHOLD=4.0 with CALLS_EXCL_FACTORIZED=2e3 -- the pT~4
# bump fix -- Y_VALS matching the plot we already made), just with a
# different LHAPDF_FILE. Running the 101 members as separate array tasks
# means the wall-clock cost is roughly ONE member's sweep, not 101x --
# they all run in parallel instead of serially like a local loop would.
#
# Output for member NNNN goes to $OUTDIR/member_NNNN/files/D0_*_LHAPDF_*.dat
# under /scratch (not /projappl -- scratch is meant for job I/O, projappl
# for the code/build). Bring the member_* directories back to this
# project's files/lhapdf/ locally (e.g. via rsync) before running
# cross_section.py's load_lhapdf_band().

set -e

# build/bin/D0 is dynamically linked against libgsl -- module load it here
# too, or every array task fails with "libgsl.so.28: cannot open shared
# object file" (the module loaded during build_roihu.sh doesn't carry over
# to a separate job's environment).
source /usr/share/lmod/lmod/init/bash
module purge
module use /appl/modulefiles
module load spack/x86_64/v2026_03/Core/gcc/15.2.0
module load gsl/2.8

LHAPDF_DIR=${LHAPDF_DIR:-data/prompt-D0-1-109}
LHAPDF_SET=${LHAPDF_SET:-prompt-D0-1-109}
OUTDIR=${OUTDIR:-/scratch/lappi/patricia/D0_lhapdf_band}

member_tag=$(printf '%04d' "$SLURM_ARRAY_TASK_ID")
member_file="$LHAPDF_DIR/${LHAPDF_SET}_${member_tag}.dat"

if [[ ! -f "$member_file" ]]; then
	echo "member $member_tag: $member_file not found -- skipping" >&2
	exit 1
fi

cd "$SLURM_SUBMIT_DIR"

echo "Task $SLURM_ARRAY_TASK_ID: member=$member_tag file=$member_file outdir=$OUTDIR/member_${member_tag}"

PT_VALS="$(seq 0.2 0.1 2.0) $(seq 2.5 0.5 12.0)" \
FRAG_TYPE=LHAPDF \
LHAPDF_FILE="$member_file" \
Y_VALS="-1.0 0.0 1.0 2.0 3.0" \
PROCESS=both \
CORES="$SLURM_CPUS_PER_TASK" \
CALLS_EXCL=1e5 CALLS_EXCL_FACTORIZED=2e3 CALLS_DIFF=1e5 \
OUTDIR="$OUTDIR/member_${member_tag}" \
./run_many_nucleus.sh

echo "Task $SLURM_ARRAY_TASK_ID done."
