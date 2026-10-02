#!/bin/bash -l
#SBATCH --job-name=D0_bk_band
#SBATCH --account=lappi
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=1000M
#SBATCH --array=0-99
#SBATCH --output=logs/bk_%A_%a.out
#SBATCH --error=logs/bk_%A_%a.err

# Usage (after build_roihu.sh and setup_bk_posterior_links.sh): sbatch run_bk_posterior_members_roihu.sbatch
# One job per BK posterior dipole sample, with the central HymnD member. Output in $OUTDIR/member_NNNN/output/.

set -e

# build/bin/D0 needs libgsl: load the module in the job too.
source /usr/share/lmod/lmod/init/bash
module purge
module use /appl/modulefiles
module load spack/x86_64/v2026_03/Core/gcc/15.2.0
module load gsl/2.8

BK_DIR=${BK_DIR:-data/Pb/bk_posterior}
OUTDIR=${OUTDIR:-/scratch/lappi/patricia/D0_bk_band}
HYMND_FILE=${HYMND_FILE:-inputs/HymnD/prompt-D0-1-109_0000.dat}

member_tag=$(printf '%04d' "$SLURM_ARRAY_TASK_ID")
member_dir="$BK_DIR/member_${member_tag}"

if [[ ! -d "$member_dir" ]]; then
	echo "member $member_tag: $member_dir not found -- run ./local_workflows/setup_bk_posterior_links.sh first" >&2
	exit 1
fi

cd "$SLURM_SUBMIT_DIR"

echo "Task $SLURM_ARRAY_TASK_ID: member=$member_tag dipole_dir=$member_dir outdir=$OUTDIR/member_${member_tag}"

PT_VALS="$(seq 0.2 0.1 2.0) $(seq 2.5 0.5 12.0)" \
DIPOLE_DIR="$member_dir" \
DIPOLE_X0=0.01 \
FRAG_TYPE=HymnD \
HYMND_FILE="$HYMND_FILE" \
Y_VALS="-1.0 0.0 1.0 2.0 3.0" \
PROCESS=both \
CORES="$SLURM_CPUS_PER_TASK" \
CALLS_EXCL=1e5 CALLS_EXCL_FACTORIZED=2e3 CALLS_DIFF=1e5 \
OUTDIR="$OUTDIR/member_${member_tag}" \
./local_workflows/run_many_nucleus.sh

echo "Task $SLURM_ARRAY_TASK_ID done."
