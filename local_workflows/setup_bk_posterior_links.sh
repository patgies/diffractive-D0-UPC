#!/bin/bash
# Usage: ./local_workflows/setup_bk_posterior_links.sh
# Makes links data/Pb/bk_posterior/member_NNNN/glauber_mve_<b> -> inputs/BK/bks_Pbtargets_1/bks/<N>/ic_208_<b>.dat.

set -e

BK_SRC=inputs/BK/bks_Pbtargets_1/bks
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
