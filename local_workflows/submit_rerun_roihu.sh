#!/bin/bash
# Submits every Pb prediction that depends on the photon flux / FF to Roihu,
# with the new effective flux (FLUX_MODEL=EFF, inputs/Gamma_AA.dat at
# sigma_NN=92 mb) and the EKO-evolved BCFY/Kniehl-Kramer grids. Same y/pT
# grids as the files currently in files/. All output goes under $BASE in
# /scratch (one subdirectory per run, see the layout at the bottom).
#
# Usage, on Roihu from the repo root:
#   srun --account=lappi --partition=small --time=00:15:00 --cpus-per-task=8 ./local_workflows/build_roihu.sh
#   STEPS=central ./local_workflows/submit_rerun_roihu.sh   # one step at a time
#   DRY_RUN=1 STEPS=central ./local_workflows/submit_rerun_roihu.sh   # only print the sbatch commands
# STEPS: any of "central scale proton xpom exclfix" (default: the first four).
# Every step writes under the same $BASE, so they can be submitted on different days.
#
# EXCL_PT_THRESHOLD=3.0 (run_many_nucleus.sh's default is 4.0): at pD0=3.5 the
# 5D exclusive integrand did not converge with 1e5 calls (negative/too-small
# values at negative y in the first central run), so the factorized integrand
# takes over from 3.0 GeV, as in run_proton_baseline.sh. "exclfix" reruns just
# the exclusive pD0=3.0, 3.5 points of that first central run this way
# (3.0 as a cross-check against the 5D result).
#
# Not rerun: D0_fixed_qp* (no photon flux, no fragmentation), ff_vs_q_scan.dat
# (HymnD only, unchanged), out/flux_scan*.csv (already current).
# Not included here: the HymnD replica band and the BK posterior band. They
# need the HymnD members 0001-0100 and the bk/ posterior dipoles on Roihu;
# once those are there, submit run_HymnD_members_roihu.sbatch /
# run_bk_posterior_members_roihu.sbatch with OUTDIR=$BASE/HymnD_band, $BASE/bk_band.

set -e

BASE=${BASE:-/scratch/lappi/patricia/rerun_EFF}
STEPS=${STEPS:-"central scale proton xpom"}
export FLUX_MODEL=EFF EXCL_PT_THRESHOLD=3.0

Y_CENTRAL="-2.0 -1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0"
PT_FINE=$(echo $(seq 0.2 0.1 2.0) $(seq 2.5 0.5 12.0))
HYMND_CENTRAL=inputs/prompt-D0-1-109/prompt-D0-1-109_0000.dat
WRAP=local_workflows/run_workflow_roihu.sbatch

if [[ ! -x build/bin/D0 || ! -x build/bin/D0_xpom ]]; then
	echo "Error: build/bin/D0 and build/bin/D0_xpom missing -- run build_roihu.sh first (see the header)." >&2
	exit 1
fi

# Run settings are exported in a subshell per job rather than passed via
# --export=VAR=..., since Y_VALS/PT_VALS contain spaces; sbatch's default
# --export=ALL hands the whole environment to the job.
submit() {
	echo "sbatch $*"
	[[ "${DRY_RUN:-0}" == "1" ]] || sbatch "$@"
}

step() { [[ " $STEPS " == *" $1 "* ]]; }

mkdir -p logs "$BASE"
echo "Output base: $BASE (FLUX_MODEL=$FLUX_MODEL, steps: $STEPS)"

# 1. Central values, one job per fragmentation function (Y_CENTRAL rapidities)
if step central; then
	for frag in BCFY KniehlKramer HymnD; do
		( export FRAG_TYPE=$frag HYMND_FILE=$HYMND_CENTRAL PROCESS=both Y_VALS="$Y_CENTRAL" PT_VALS="$PT_FINE" \
		         CALLS_EXCL=1e5 CALLS_EXCL_FACTORIZED=2e3 CALLS_DIFF=1e5 OUTDIR=$BASE/central
		  submit --job-name="D0_central_${frag}" $WRAP local_workflows/run_many_nucleus.sh )
	done
fi

# 1b. Exclusive pD0=3.0, 3.5 with the factorized integrand, to patch the first
#     central run (made with EXCL_PT_THRESHOLD=4.0); see the header
#     EXCLFIX_PT / EXCLFIX_DIR pick other points to check, e.g.
#     STEPS=exclfix EXCLFIX_PT=2.5 EXCLFIX_DIR=central_exclfix_pt25
if step exclfix; then
	EXCLFIX_PT=${EXCLFIX_PT:-"3.0 3.5"}
	EXCLFIX_DIR=${EXCLFIX_DIR:-central_exclfix}
	thr=$(echo $EXCLFIX_PT | awk '{m=$1; for (i=2;i<=NF;i++) if ($i<m) m=$i; print m}')   # factorized at every point asked for
	for frag in BCFY KniehlKramer HymnD; do
		( export FRAG_TYPE=$frag HYMND_FILE=$HYMND_CENTRAL PROCESS=exclusive Y_VALS="$Y_CENTRAL" PT_VALS="$EXCLFIX_PT" \
		         EXCL_PT_THRESHOLD=$thr CALLS_EXCL_FACTORIZED=2e3 OUTDIR=$BASE/$EXCLFIX_DIR
		  submit --job-name="D0_exclfix_${frag}" --time=01:00:00 --cpus-per-task=32 $WRAP local_workflows/run_many_nucleus.sh )
	done
fi

# 2. Fragmentation-scale variation (Q = 0.5 and 2 m_T, Y_CENTRAL rapidities)
#    for all three FFs -- BCFY and Kniehl-Kramer are DGLAP-evolved too (EKO
#    grids, Q = 1.5-50 GeV), so they get the same envelope as HymnD
if step scale; then
	for frag in HymnD BCFY KniehlKramer; do
		for factor in 0.5 2.0; do
			( export FRAG_TYPE=$frag HYMND_FILE=$HYMND_CENTRAL SCALE_FACTOR=$factor PROCESS=both Y_VALS="$Y_CENTRAL" PT_VALS="$PT_FINE" \
			         CALLS_EXCL=1e5 CALLS_EXCL_FACTORIZED=2e3 CALLS_DIFF=1e5 OUTDIR=$BASE/${frag}_scale/factor_${factor}
			  submit --job-name="D0_scale_${frag}_${factor}" $WRAP local_workflows/run_many_nucleus.sh )
		done
	done
fi

# 3. Proton baseline for R_pA (same FLUX_MODEL as the Pb runs; no b loop, cheap)
if step proton; then
	( export FRAG_TYPE=HymnD HYMND_FILE=$HYMND_CENTRAL OUTDIR=$BASE/proton_baseline
	  submit --job-name=D0_proton_baseline --time=01:00:00 --cpus-per-task=20 $WRAP local_workflows/run_proton_baseline.sh )
fi

# 4. x_po spectrum: BCFY at pD0=2 GeV (y=0-3), what python/xpom_spectrum.py
#    plots by default. (The old Kniehl-Kramer x_po set is not rerun; for it:
#    FRAG_TYPE=KniehlKramer with run_many_xpom.sh's default grid, ~64 cores x 8 h.)
if step xpom; then
	( export FRAG_TYPE=BCFY Y_VALS="0.0 1.0 2.0 3.0" PT_VALS="2.0" OUTDIR=$BASE/xpom
	  submit --job-name=D0_xpom_BCFY --time=02:00:00 --cpus-per-task=32 $WRAP local_workflows/run_many_xpom.sh )
fi

cat <<EOF

Submitted ($STEPS). When everything has finished, $BASE contains:
  central/files/D0_{exclusive,diffractive}_{BCFY,KniehlKramer,HymnD}_An0n_Pb_y*.dat  -> files/
  HymnD_scale/factor_{0.5,2.0}/files/                                              -> files/HymnD_scale/factor_*/files/
  {BCFY,KniehlKramer}_scale/factor_{0.5,2.0}/files/                                -> files/{BCFY,KniehlKramer}_scale/factor_*/files/
  proton_baseline/files/D0_proton_baseline_*                                       -> files/
  xpom/files/D0_*_xpom_BCFY_*                                                      -> files/
then rebuild combined_d0_* with python/combine_nucleus.py.
EOF
