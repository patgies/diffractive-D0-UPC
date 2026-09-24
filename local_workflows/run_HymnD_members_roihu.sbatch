#!/bin/bash -l
#SBATCH --job-name=D0_HymnD_band
#SBATCH --account=lappi
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=1000M
#SBATCH --array=0-100
#SBATCH --output=logs/HymnD_%A_%a.out
#SBATCH --error=logs/HymnD_%A_%a.err

# Usage (run build_roihu.sh first -- this script does NOT build):
#   sbatch run_HymnD_members_roihu.sbatch
#   sbatch --array=0-9 run_HymnD_members_roihu.sbatch          # quick subset test
#   sbatch --export=ALL,HYMND_DIR=/other/path run_HymnD_members_roihu.sbatch
#   # one-off central-value run at extra rapidities (member 0000 only,
#   # e.g. filling in more y points for a contour plot -- not a replica
#   # sweep, so --array=0-0 and a different OUTDIR so it doesn't land
#   # inside the replica-band output):
#   sbatch --array=0-0 --export=ALL,Y_VALS="-3.0 -2.5 -2.0 -1.5 -0.5 0.5 1.5 2.5",OUTDIR=/scratch/lappi/patricia/D0_contour run_HymnD_members_roihu.sbatch
#   # one-off run with a different fragmentation scheme entirely (BCFY,
#   # KniehlKramer -- no HYMND_FILE/member lookup happens for these, so
#   # --array=0-0 is enough regardless of member count):
#   sbatch --array=0-0 --export=ALL,FRAG_TYPE=BCFY,Y_VALS="-3.0 -2.5 -2.0 -1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0 2.5 3.0",OUTDIR=/scratch/lappi/patricia/D0_bcfy run_HymnD_members_roihu.sbatch
#
# One array task per HymnD replica member (0-100: member 0 is the central
# value, 1-100 are the replicas -- see inputs/prompt-D0-1-109/*.info,
# ErrorType=replicas). Each task runs the exact same diffractive+exclusive
# sweep as the local member-0000 run (fine pT grid below 2 GeV, coarser
# above, EXCL_PT_THRESHOLD=4.0 with CALLS_EXCL_FACTORIZED=2e3 -- the pT~4
# bump fix -- Y_VALS matching the plot we already made), just with a
# different HYMND_FILE. Running the 101 members as separate array tasks
# means the wall-clock cost is roughly ONE member's sweep, not 101x --
# they all run in parallel instead of serially like a local loop would.
#
# Output for member NNNN goes to $OUTDIR/member_NNNN/files/D0_*_HymnD_*.dat
# under /scratch (not /projappl -- scratch is meant for job I/O, projappl
# for the code/build). Bring the member_* directories back to this
# project's files/HymnD/ locally (e.g. via rsync) before running
# cross_section.py's load_hymnd_replica_band().

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

HYMND_DIR=${HYMND_DIR:-inputs/prompt-D0-1-109}
HYMND_SET=${HYMND_SET:-prompt-D0-1-109}
OUTDIR=${OUTDIR:-/scratch/lappi/patricia/D0_HymnD_band}
Y_VALS=${Y_VALS:-"-1.0 0.0 1.0 2.0 3.0"}
FRAG_TYPE=${FRAG_TYPE:-HymnD}

member_tag=$(printf '%04d' "$SLURM_ARRAY_TASK_ID")
member_file=""
if [[ "$FRAG_TYPE" == "HymnD" ]]; then
	member_file="$HYMND_DIR/${HYMND_SET}_${member_tag}.dat"
	if [[ ! -f "$member_file" ]]; then
		echo "member $member_tag: $member_file not found -- skipping" >&2
		exit 1
	fi
fi

cd "$SLURM_SUBMIT_DIR"

echo "Task $SLURM_ARRAY_TASK_ID: frag=$FRAG_TYPE member=$member_tag file=$member_file outdir=$OUTDIR/member_${member_tag}"

PT_VALS="$(seq 0.2 0.1 2.0) $(seq 2.5 0.5 12.0)" \
FRAG_TYPE="$FRAG_TYPE" \
HYMND_FILE="$member_file" \
Y_VALS="$Y_VALS" \
PROCESS=both \
CORES="$SLURM_CPUS_PER_TASK" \
CALLS_EXCL=1e5 CALLS_EXCL_FACTORIZED=2e3 CALLS_DIFF=1e5 \
OUTDIR="$OUTDIR/member_${member_tag}" \
./local_workflows/run_many_nucleus.sh

echo "Task $SLURM_ARRAY_TASK_ID done."
