#!/bin/bash
ensure_exists() {
	for f in "$@"; do
		if ! [ -e "$f" ]; then
			echo missing "$f"
			result=1
		fi
	done
}

ensure_count() {
	if [ $(ls $1 | wc -l) -lt "$2" ]; then
		echo missing "$1" got $(ls $1 | wc -l) instead of "$2"
		result=1
	fi
}

n_snap=91

result=0
ensure_count "subfind/fof_subhalo_tab_*" $n_snap
ensure_exists sublink/tree.hdf5 sublink/tree_extended.hdf5
ensure_count "sublink/offsets/offsets_*.hdf5" $n_snap
ensure_exists rockstar/trees/tree_0_0_0.dat rockstar/trees/locations.dat rockstar/trees/forests.list
ensure_count rockstar/hlists/ $(($n_snap - 20)) # This isn't an exact count since it'll skip snapshots with not enough halos
ensure_count CMD/2D_maps 12
ensure_count CMD/3D_grids $((3 * 5 * 12))
ensure_count Pk/ $((($n_snap + 1) * 5))
ensure_exists DisPerSE/massgrid-disperse_G-256_S-02_c500-output_file_090.hdf5

if [ "$result" -eq 0 ]; then
	true
else
	false
fi
