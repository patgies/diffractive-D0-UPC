#!/bin/bash
# Runs run_nucleus.sh for each HymnD member or BK sample, in output/<CHANNEL>/<set>_band/member_NNNN/.
# Usage: MEMBER_SET=HymnD (or bk) MEMBERS="0 1 2" ./run_scripts/run_members.sh

set -e
cd "$(dirname "$0")/.."
source run_scripts/config.sh

require_file() { [[ -f "$1" ]] || { echo "Error: $1 does not exist. $2" >&2; exit 1; }; }
require_dir()  { [[ -d "$1" ]] || { echo "Error: $1 does not exist. $2" >&2; exit 1; }; }

MEMBER_SET=${MEMBER_SET:-HymnD}
HYMND_DIR=${HYMND_DIR:-input/HymnD}
HYMND_SET=${HYMND_SET:-prompt-D0-1-109}
BK_DIR=${BK_DIR:-data/Pb/bk_posterior}

case "$MEMBER_SET" in
	HymnD)
		OUTBASE=${OUTBASE:-$OUTPUT_ROOT/$channel_tag/HymnD_band}
		all_members=$(ls "$HYMND_DIR"/${HYMND_SET}_[0-9][0-9][0-9][0-9].dat 2>/dev/null | sed -E 's/.*_([0-9]{4})\.dat/\1/')
		;;
	bk)
		OUTBASE=${OUTBASE:-$OUTPUT_ROOT/$channel_tag/bk_band}
		all_members=$(ls -d "$BK_DIR"/member_[0-9][0-9][0-9][0-9] 2>/dev/null | sed -E 's/.*member_([0-9]{4})/\1/')
		;;
	*)
		echo "Error: MEMBER_SET must be HymnD or bk, not '$MEMBER_SET'." >&2
		exit 1
		;;
esac

MEMBERS=${MEMBERS:-$all_members}
if [[ -z "$MEMBERS" ]]; then
	echo "Error: no $MEMBER_SET members found (see the header of this script)." >&2
	exit 1
fi

for m in $MEMBERS; do
	member_tag=$(printf '%04d' "$((10#$m))")
	echo "=== $MEMBER_SET member $member_tag ($(date)) ==="
	if [[ "$MEMBER_SET" == "HymnD" ]]; then
		member_file="$HYMND_DIR/${HYMND_SET}_${member_tag}.dat"
		require_file "$member_file"
		FRAG_TYPE=HymnD HYMND_FILE="$member_file" OUTDIR="$OUTBASE/member_${member_tag}" \
			./run_scripts/run_nucleus.sh
	else
		member_dir="$BK_DIR/member_${member_tag}"
		require_dir "$member_dir" "Run ./run_scripts/setup_bk_posterior_links.sh first."
		DIPOLE_DIR="$member_dir" DIPOLE_X0=0.01 OUTDIR="$OUTBASE/member_${member_tag}" \
			./run_scripts/run_nucleus.sh
	fi
done

echo "Finished at $(date)"
