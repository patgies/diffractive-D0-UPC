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

# Usage (run build_roihu.sh and ./local_workflows/setup_bk_posterior_links.sh first --
# this script does NOT build or set up symlinks):
#   sbatch run_bk_posterior_members_roihu.sbatch
#   sbatch --array=0-9 run_bk_posterior_members_roihu.sbatch   # quick subset test
#
# Counterpart to run_HymnD_members_roihu.sbatch, but for the *other*
# uncertainty source: instead of varying the fragmentation function while
# holding the dipole amplitude fixed, this holds the fragmentation scheme
# fixed (FRAG_TYPE=HymnD, central member 0000 -- matching the plotted
# central HymnD curve) and varies the dipole amplitude across 100 posterior
# samples of the BK initial-condition fit ("LOmvefit": Q_{s,0}^2, e_c, C^2,
# sigma0/2 -- see bk/posteriorsamples_100_LOmvefit.dat). This isolates the
# BK-initial-condition uncertainty as a second, independent band -- meant to
# be combined with the HymnD-replica band in quadrature, not cross-produced
# with it (that would be 100x101 runs for no real precision gain).
#
# data/Pb/bk_posterior/member_<NNNN>/glauber_mve_<b> are symlinks into
# bk/bks_Pbtargets_1/bks/<NNNN>/ic_208_<b>.dat (see
# setup_bk_posterior_links.sh). Those files leave x0 out of their header,
# which DataFile reads as "invalid x0" and resets to 0 -- DIPOLE_X0=0.01
# below works around that via src/main.cpp's SetX0() override (0.01 matches
# both data/Pb/mve/'s header and the standard HERA-fit convention this
# family of Bayesian MVe fits uses).
#
# Output for member NNNN goes to $OUTDIR/member_NNNN/output/D0_*_HymnD_*.dat
# under /scratch. Bring the member_* directories back to this project's
# output/bk_posterior/ locally (e.g. via rsync) before running
# D0.py's load_bk_band().

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

BK_DIR=${BK_DIR:-data/Pb/bk_posterior}
OUTDIR=${OUTDIR:-/scratch/lappi/patricia/D0_bk_band}
HYMND_FILE=${HYMND_FILE:-inputs/prompt-D0-1-109/prompt-D0-1-109_0000.dat}

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
