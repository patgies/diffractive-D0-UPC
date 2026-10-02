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

# Usage (after build_roihu.sh): sbatch run_HymnD_members_roihu.sbatch   [--array=0-9 for a subset]
# One job per HymnD member (0 = central, 1-100 = replicas). Output in $OUTDIR/member_NNNN/output/.

set -e

# build/bin/D0 needs libgsl: load the module in the job too.
source /usr/share/lmod/lmod/init/bash
module purge
module use /appl/modulefiles
module load spack/x86_64/v2026_03/Core/gcc/15.2.0
module load gsl/2.8

HYMND_DIR=${HYMND_DIR:-inputs/HymnD}
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
