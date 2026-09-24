#!/bin/bash
# Usage: ./local_workflows/setup_bk_posterior_links.sh
#
# The BK-initial-condition posterior sample dipole files live at
# bk/bks_Pbtargets_1/bks/<N>/ic_208_<b>.dat (N = 0..99, b = 0,2,...,34 -- see
# bk/posteriorsamples_100_LOmvefit.dat, one row per N). run_many_nucleus.sh
# (and the underlying D0 binary) expect a DIPOLE_DIR containing files named
# glauber_mve_<b>, matching data/Pb/mve/'s convention -- so this creates
# relative symlinks data/Pb/bk_posterior/member_NNNN/glauber_mve_<b> ->
# ../../../../bk/bks_Pbtargets_1/bks/<N>/ic_208_<b>.dat for every member,
# rather than duplicating the 359M of underlying bk/ data or changing
# run_many_nucleus.sh's glob pattern. Safe to re-run (skips existing links).
#
# These files also leave x0 out of their header (DataFile reads that as
# "invalid x0" and resets it to 0) -- run_bk_posterior_members_roihu.sbatch
# passes DIPOLE_X0=0.01 to work around that, see src/main.cpp.

set -e

BK_SRC=bk/bks_Pbtargets_1/bks
BK_DIR=data/Pb/bk_posterior

for member_dir in "$BK_SRC"/*; do
	[[ -d "$member_dir" ]] || continue
	n=$(basename "$member_dir")
	member_tag=$(printf '%04d' "$n")
	dest="$BK_DIR/member_${member_tag}"
	mkdir -p "$dest"
	for f in "$member_dir"/ic_208_*.dat; do
		b=$(basename "$f" | sed -E 's/ic_208_([0-9]+)\.dat/\1/')
		ln -srf "$f" "$dest/glauber_mve_${b}"
	done
done

n_members=$(ls -d "$BK_DIR"/member_* | wc -l)
echo "Done: $n_members member directories under $BK_DIR/"
