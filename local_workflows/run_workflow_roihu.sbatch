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

# Runs one workflow script as a Roihu job. Settings are given as environment variables:
#   sbatch --job-name=<name> --export=ALL,FRAG_TYPE=BCFY,OUTDIR=<dir> local_workflows/run_workflow_roihu.sbatch <script>

set -e

if [[ $# -lt 1 ]]; then
	echo "Usage: sbatch [options] run_workflow_roihu.sbatch <workflow_script>" >&2
	exit 1
fi

# build/bin/D0 needs libgsl: load the module in the job too
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
