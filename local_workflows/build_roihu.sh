#!/bin/bash
# Build inside a job, on the same partition where you run (not on the login node):
#   srun --account=lappi --partition=small --time=00:15:00 --cpus-per-task=8 ./local_workflows/build_roihu.sh

set -e

source /usr/share/lmod/lmod/init/bash
module purge
module use /appl/modulefiles
module load spack/x86_64/v2026_03/Core/gcc/15.2.0
module load gsl/2.8

mkdir -p build
cmake -S . -B build
cmake --build build -j"$(nproc)" --target D0 D0_xpom

echo "Build OK: build/bin/D0, build/bin/D0_xpom"
