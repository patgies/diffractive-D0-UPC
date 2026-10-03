#!/bin/bash
# Default settings of all workflows. Change any of them on the command line,
# e.g. CHANNEL=0n0n ./run_scripts/run_nucleus.sh

# Collision system and photon flux
: "${NUCLEUS:=Pb}"                 # target nucleus
: "${CHANNEL:=An0n}"               # An0n | Xn0n | 0n0n | AnAn | PL(AnAn)
: "${FLUX_MODEL:=EFF}"             # EFF | STARLIGHT | PL | WS
: "${TARGET:=AA}"                  # AA (Pb+Pb, 5.36 TeV) | pA (p+Pb, 8.16 TeV)

# Fragmentation
: "${FRAG_TYPE:=HymnD}"            # BCFY | KniehlKramer | HymnD
: "${HYMND_FILE:=input/HymnD/prompt-D0-1-109_0000.dat}"   # central HymnD member
: "${SCALE_FACTOR:=1.0}"           # fragmentation scale Q = SCALE_FACTOR * m_T
: "${SCALE_FACTORS:=0.5 2.0}"      # scale factors of run_scale_variation.sh

# Kinematic grid
: "${Y_VALS:=-2.0 -1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0}"
: "${PT_VALS:=$(echo $(seq 0.2 0.1 2.0) $(seq 2.5 0.5 12.0))}"
: "${PROCESS:=both}"               # exclusive | diffractive | both

# VEGAS calls per point. Above EXCL_PT_THRESHOLD the exclusive part uses the
# faster 3D integral, with CALLS_EXCL_FACTORIZED calls.
: "${CALLS:=1e5}"
: "${CALLS_EXCL:=$CALLS}"
: "${CALLS_DIFF:=$CALLS}"
: "${CALLS_EXCL_FACTORIZED:=2e3}"
: "${EXCL_PT_THRESHOLD:=3.0}"

# x_po grid of run_xpom.sh (log-spaced)
: "${XPO_MIN:=1e-6}"
: "${XPO_MAX:=0.1}"
: "${XPO_N:=25}"

# Folder for all results
: "${OUTPUT_ROOT:=output}"

# Parallel processes (in a cluster job, set it to the number of CPUs)
: "${CORES:=$(( $(nproc) / 2 ))}"

export OUTPUT_ROOT NUCLEUS CHANNEL FLUX_MODEL TARGET FRAG_TYPE HYMND_FILE SCALE_FACTOR PROCESS
export CALLS CALLS_EXCL CALLS_DIFF CALLS_EXCL_FACTORIZED EXCL_PT_THRESHOLD

# Tags used in folder and file names
channel_tag=$(echo "$CHANNEL" | tr -d '() ')
# empty for EFF, "_STARLIGHT" etc. for the other fluxes
flux_tag=$([[ "$FLUX_MODEL" == "EFF" ]] && echo "" || echo "_${FLUX_MODEL}")
