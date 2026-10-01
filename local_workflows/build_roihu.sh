#!/bin/bash
# CMakeLists.txt sets -march=native, and unlike Puhti (where login and
# compute nodes were the same Cascade Lake generation), Roihu's login node
# (AMD EPYC 9335, 32-core) does NOT look like the same chip as the "small"
# partition's compute nodes (384 logical CPUs/node = a much bigger part,
# possibly a newer generation) -- so -march=native compiled here could pick
# instructions the compute nodes don't have, or vice versa. Don't build on
# the login node; build inside a quick job on the same partition you'll run
# on instead:
#
#   srun --account=lappi --partition=small --time=00:15:00 --cpus-per-task=8 ./local_workflows/build_roihu.sh
#
# Roihu's system gcc/cmake (in /usr/bin) are fine to build with, but GSL
# isn't installed system-wide -- it comes from the Spack module tree under
# /appl/modulefiles, which is hierarchical: a Core compiler module has to be
# loaded before "gsl" becomes available at all.
#
# NOTE: exact module names/versions drift over time -- if `module load`
# below fails, run `module use /appl/modulefiles && module spider gsl` to
# see what's currently installed and adjust.

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
