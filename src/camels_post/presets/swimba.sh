# This function should take a number and return the path to the a snapshot of that number.
# The snapshot should be in a gadget-hdf5 compatible format.
# Note that many aspects of this pipeline expect snapshot names to all be of the form `some/path/to/snapshotFileName_<snapNumber>.hdf5`
# If this is not the case, you may have to alter things further in the code to make it work.
# I can probably set this up such that it symlinks files to be the correct format, bypassing this requirement.
function get-gadget-snapshot() {
	printf "/scratch/${SIM_ROOT}/snap_%03d.hdf5" "$1"
}

function get-gadget-ic() {
	echo "/scratch/${SIM_ROOT}/ic.hdf5"
}

# This step is necessary for SWIMBA because many of the steps don't work natively with swift-style snapshots.
# If your snapshots are already gadget-hdf5 compatible, then this part can be omitted. Also you should turn off
# CLEAN_GADGET_SNAPS and CONVERT_SNAPSHOTS
function convert-snaps() {
	mkdir -p "/scratch/${SIM_ROOT}/"

	# I'll grant this is a pretty roundabout way to do this.
	# For now though this is effectively a for loop going through every snapshot.
	# It's just weird because I'm using xargs to parallelize the process.
	# Also, swift2gadget is included in my Swimba Pipeline, not my postprocessing package.
	seq 0 90 | xargs -P"$cpus" -I{} bash -c '
        SOURCE_SNAP=$(printf "snaps/snapshot_%04d.hdf5" {})
        DEST_SNAP=$(printf "/scratch/${SIM_ROOT}/snap_%03d.hdf5" {})

        # Does the gadget snapshot already exist?
        if [ -e "${DEST_SNAP}" ]; then
            # Skip it!
            exit
        fi

        swift2gadget "${SOURCE_SNAP}" "${DEST_SNAP}"
    '

	if [ ! -e "/scratch/${SIM_ROOT}/ic.hdf5" ]; then
		combine-IC ICs/ics "/scratch/${SIM_ROOT}/ic.hdf5"
	fi
}

# Any additional cleanup which needs to happen.
# For example, in SWIMBA I put the converted snapshots in node-local scratch storage.
# Slurm on rusty will clean this up itself, but it's better to clean it up myself.
function cleanup() {
	rm -r "/scratch/${SIM_ROOT}/"
}

# Get softening in kpc / h
# This is necessary for subfind and rockstar, since subhalo finding depends on softening lengths.
# $1 is a particle type
# $2 is either "comoving" or "physical"
# If the softening length does not change across particle type or the max physical softening length is equal to the comoving softening length,
# then this can just ignore the arguments and return that one value.
# This implementation just extracts the corresponding values from the SWIFT parameter file.
# Note the fact / h!!!
# For SWIFT at least, I've got to also extract h from the parameter file and perform that conversion.
# At the very least, it should support particle types 0-4 since those are expected to be in the Subfind parameter file.
function get-softening() {
	h=$(sed -ne 's/^\s*h: //p' params.yml)
	l=$(
		if [ "$2" = "comoving" ]; then
			case $1 in
			1 | 2)
				sed -ne 's/^\s*comoving_DM_softening: //p' params.yml
				;;
			0 | 3 | 4 | 5)
				sed -ne 's/^\s*comoving_baryon_softening: //p' params.yml
				;;
			esac
		elif [ "$2" = "physical" ]; then
			case $1 in
			1 | 2)
				sed -ne 's/^\s*max_physical_DM_softening: //p' params.yml
				;;
			0 | 3 | 4 | 5)
				sed -ne 's/^\s*max_physical_baryon_softening: //p' params.yml
				;;
			esac
		fi
	)
	bc <<<"$l * $h * 1000"
}

# ALL_SNAPS should be an array containing every snapshot number we're outputting.
# I've set it up like this to allow for non-sequential snapshot numbers.
mapfile -t ALL_SNAPS < <(seq 0 90)

# The _OUTPUT variables are for determining where the ouptuts for the various steps should go.
# Note that some steps depend on previous steps, so changing one output can break later steps.
SUBFIND_OUTPUT=./subfind/
SUBLINK_OUTPUT=./sublink/
SUBLINK_GAL_OUTPUT=./sublink_gal/
ROCKSTAR_OUTPUT=./rockstar/
CMD_OUTPUT=./CMD/
PK_OUTPUT=./Pk/
DISPERSE_OUTPUT=./DisPerSE/

# The _ROOT variables are for where the various codes live.
# There are four version of Sublink:
# - 4 byte particle IDs, single-file snapshots
# - 4 byte particle IDs, multi-file snapshots
# - 8 byte particle IDs, single-file snapshots
# - 8 byte particle IDs, multi-file snapshots
# This is the first of those versions.
# You should change SUBLINK_ROOT to point to the corresponding version you need.
SUBLINK_ROOT=/mnt/home/ssutherland/codes/SubLink/SubLink_SHORTIDS/
AREPO_ROOT=/mnt/home/ssutherland/codes/Arepo_subfind_v2/
ROCKSTAR_ROOT=/mnt/ceph/users/camels/Codes/rockstar-galaxies/
CONSISTENT_TREE_ROOT=/mnt/ceph/users/camels/Codes/consistent-trees/
DISPERSE_ROOT=/mnt/home/ssutherland/codes/DisPerSE/bin/

# In the CAMELS Multifield Dataset, we make 3D grids for z = 0, 0.5, 1, 1.5, and 2.
# This should be an array with the snapshot numbers corresponding with these redshifts.
CMD_SNAPSHOTS=(90 73 61 52 44)

# When the -f option is passed, this script will delete all the existing data products.
# If this is "yes", then doing so will also delete the gadget snapshots.
# IF THE GADGET SNAPSHOTS RETURNED BY get-gadget-snapshot ARE THE ORIGINAL DATA, MAKE SURE THIS IS UNSET.
# OTHERWISE IT WILL DELETE THE ORIGINAL SNAPSHOTS.
CLEAN_GADGET_SNAPS=yes
# Should convert-snapshots be run before the rest of the pipeline does?
CONVERT_SNAPSHOTS=yes
# To get CAMELS parameters in the header for the DisPerSE catalogs, this code uses the cosmoastroseed parameter files present in the camels directories.
# If this file is not present, or you simply don't care about having the parameter values in the header, then you can turn this off.
WITH_COSMOASTROSEED=yes
