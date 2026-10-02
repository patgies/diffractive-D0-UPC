#!/bin/bash
# Submits the Pb runs to Roihu. Output goes to $BASE. STEPS: central scale proton xpom exclfix.
# Usage (repo root, after build_roihu.sh): STEPS=central ./local_workflows/submit_rerun_roihu.sh   [DRY_RUN=1]

set -e

BASE=${BASE:-/scratch/lappi/patricia/rerun_EFF}
STEPS=${STEPS:-"central scale proton xpom"}
# EXCL_PT_THRESHOLD=3.0: the 5D exclusive integral was not stable at pD0=3.5 with 1e5 calls
export FLUX_MODEL=EFF EXCL_PT_THRESHOLD=3.0

Y_CENTRAL="-2.0 -1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0"
PT_FINE=$(echo $(seq 0.2 0.1 2.0) $(seq 2.5 0.5 12.0))
HYMND_CENTRAL=inputs/HymnD/prompt-D0-1-109_0000.dat
WRAP=local_workflows/run_workflow_roihu.sbatch

if [[ ! -x build/bin/D0 || ! -x build/bin/D0_xpom ]]; then
	echo "Error: build/bin/D0 and build/bin/D0_xpom missing -- run build_roihu.sh first (see the header)." >&2
	exit 1
fi

# Settings are exported for each job in a subshell (Y_VALS and PT_VALS have spaces).
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

# 1b. Exclusive at EXCLFIX_PT (default 3.0, 3.5) with the 3D version. Output in EXCLFIX_DIR.
if step exclfix; then
	EXCLFIX_PT=${EXCLFIX_PT:-"3.0 3.5"}
	EXCLFIX_DIR=${EXCLFIX_DIR:-central_exclfix}
	thr=$(echo $EXCLFIX_PT | awk '{m=$1; for (i=2;i<=NF;i++) if ($i<m) m=$i; print m}')   # use the 3D version at every point asked for
	for frag in BCFY KniehlKramer HymnD; do
		( export FRAG_TYPE=$frag HYMND_FILE=$HYMND_CENTRAL PROCESS=exclusive Y_VALS="$Y_CENTRAL" PT_VALS="$EXCLFIX_PT" \
		         EXCL_PT_THRESHOLD=$thr CALLS_EXCL_FACTORIZED=2e3 OUTDIR=$BASE/$EXCLFIX_DIR
		  submit --job-name="D0_exclfix_${frag}" --time=01:00:00 --cpus-per-task=32 $WRAP local_workflows/run_many_nucleus.sh )
	done
fi

# 2. Fragmentation-scale variation (Q = 0.5 and 2 m_T) for all three FFs
if step scale; then
	for frag in HymnD BCFY KniehlKramer; do
		for factor in 0.5 2.0; do
			( export FRAG_TYPE=$frag HYMND_FILE=$HYMND_CENTRAL SCALE_FACTOR=$factor PROCESS=both Y_VALS="$Y_CENTRAL" PT_VALS="$PT_FINE" \
			         CALLS_EXCL=1e5 CALLS_EXCL_FACTORIZED=2e3 CALLS_DIFF=1e5 OUTDIR=$BASE/${frag}_scale/factor_${factor}
			  submit --job-name="D0_scale_${frag}_${factor}" $WRAP local_workflows/run_many_nucleus.sh )
		done
	done
fi

# 3. Proton baseline for R_pA (same FLUX_MODEL as the Pb runs, no b loop, fast)
if step proton; then
	( export FRAG_TYPE=HymnD HYMND_FILE=$HYMND_CENTRAL OUTDIR=$BASE/proton_baseline
	  submit --job-name=D0_proton_baseline --time=01:00:00 --cpus-per-task=20 $WRAP local_workflows/run_proton_baseline.sh )
fi

# 4. x_po spectrum: BCFY at pD0 = 2 GeV, y = 0-3 (what python/xpom.py plots)
if step xpom; then
	( export FRAG_TYPE=BCFY Y_VALS="0.0 1.0 2.0 3.0" PT_VALS="2.0" OUTDIR=$BASE/xpom
	  submit --job-name=D0_xpom_BCFY --time=02:00:00 --cpus-per-task=32 $WRAP local_workflows/run_many_xpom.sh )
fi

cat <<EOF

Submitted ($STEPS). When everything has finished, $BASE contains:
  central/output/D0_{exclusive,diffractive}_{BCFY,KniehlKramer,HymnD}_An0n_Pb_y*.dat  -> output/
  HymnD_scale/factor_{0.5,2.0}/output/                                              -> output/HymnD_scale/factor_*/output/
  {BCFY,KniehlKramer}_scale/factor_{0.5,2.0}/output/                                -> output/{BCFY,KniehlKramer}_scale/factor_*/output/
  proton_baseline/output/D0_proton_baseline_*                                       -> output/
  xpom/output/D0_*_xpom_BCFY_*                                                      -> output/
then rebuild combined_d0_* with python/combine_nucleus.py.
EOF
