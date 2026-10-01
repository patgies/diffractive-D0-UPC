#!/bin/bash -l
#SBATCH --job-name=D0_workflow
#SBATCH --account=lappi
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=1000M
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Generic Roihu wrapper: runs ONE local workflow script (run_many_nucleus.sh,
# run_proton_baseline.sh, run_many_xpom.sh, ...) as a single job, with the
# module environment loaded and CORES set to the job's CPU count. All the
# settings (FRAG_TYPE, Y_VALS, OUTDIR, FLUX_MODEL, ...) are passed as
# environment variables, exactly as when running the script locally:
#
#   sbatch --job-name=D0_central_BCFY \
#          --export=ALL,FRAG_TYPE=BCFY,PROCESS=both,OUTDIR=/scratch/lappi/patricia/xyz \
#          local_workflows/run_workflow_roihu.sbatch local_workflows/run_many_nucleus.sh
#
# --time/--cpus-per-task given on the sbatch command line override the
# defaults above. Run from the repo root, after build_roihu.sh (this does
# NOT build; SKIP_BUILD=1 below stops run_proton_baseline.sh/run_many_xpom.sh
# from running cmake in parallel jobs on the same build/ directory).
# See submit_rerun_roihu.sh for the full set of production runs.

set -e

if [[ $# -lt 1 ]]; then
	echo "Usage: sbatch [options] run_workflow_roihu.sbatch <workflow_script>" >&2
	exit 1
fi

# build/bin/D0 is dynamically linked against libgsl -- load it in the job too
source /usr/share/lmod/lmod/init/bash
module purge
module use /appl/modulefiles
module load spack/x86_64/v2026_03/Core/gcc/15.2.0
module load gsl/2.8

cd "$SLURM_SUBMIT_DIR"
export CORES="$SLURM_CPUS_PER_TASK"
export SKIP_BUILD=1

echo "Job $SLURM_JOB_ID ($SLURM_JOB_NAME): $1, OUTDIR=${OUTDIR:-.}, FLUX_MODEL=${FLUX_MODEL:-EFF}, $(date)"
bash "$1"
echo "Job $SLURM_JOB_ID done, $(date)"
