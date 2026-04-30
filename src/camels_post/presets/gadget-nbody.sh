# This function should take a number and return the path to the a snapshot of that number.
# The snapshot should be in a gadget-hdf5 compatible format.
# Note that many aspects of this pipeline expect snapshot names to all be of the form `some/path/to/snapshotFileName_<snapNumber>.hdf5`
# If this is not the case, you may have to alter things further in the code to make it work.
# I can probably set this up such that it symlinks files to be the correct format, bypassing this requirement.
function get-gadget-snapshot() {
	printf "./snap_%03d.hdf5" "$1"
}

function get-gadget-ic() {
	echo "./ICs/ic_combined.hdf5"
}

# This step is necessary for SWIMBA because many of the steps don't work natively with swift-style snapshots.
# If your snapshots are already gadget-hdf5 compatible, then this part can be omitted. Also you should turn off
# CLEAN_GADGET_SNAPS and CONVERT_SNAPSHOTS
function convert-snaps() {
	if [ ! -e "./ICs/ic_combined.hdf5" ]; then
		combine-IC "ICs/ics" "ICs/ic_combined.hdf5"
	fi
}

# Any additional cleanup which needs to happen.
# For example, in SWIMBA I put the converted snapshots in node-local scratch storage.
# Slurm on rusty will clean this up itself, but it's better to clean it up myself.
function cleanup() {
	:
}

# Get softening in Mpc / h
# This is necessary for subfind and rockstar, since subhalo finding depends on softening lengths.
# $1 is a particle type
# $2 is either "comoving" or "physical"
# If the softening length does not change across particle type or the max physical softening length is equal to the comoving softening length,
# then this can just ignore the arguments and return that one value.
# This implementation just extracts the corresponding values from the SWIFT parameter file.
# Note the fact / h!!!
# For SWIFT at least, I've got to also extract h from the parameter file and perform that conversion.
# At the very least, it should support particle types 0-4 since those are expected to be in the Subfind parameter file.
# TODO: CHECK SOFTENING UNITS FOR SUBFIND
function get-softening() {
	case $1 in
	0)
		type=Gas
		;;
	1)
		type=Halo
		;;
	2)
		type=Disk
		;;
	3)
		type=Bulge
		;;
	4)
		type=Stars
		;;
	esac
	if [ "$2" = "comoving" ]; then
		suffix="MaxPhys"
	else
		suffix=""
	fi
	param="Softening${type}${suffix}"
	# -w needed to prevent SofteningGasMaxPhys from showing up in SofteningGas
	grep -w "$param" parameters-usedvalues | awk '{print($2)}'
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

# Particle types which are baryons
BARYON_PTYPES=(0 4 5)
# Particle types which are dark matter
DM_PTYPES=(1)
# In the CAMELS Multifield Dataset, we make 3D grids for z = 0, 0.5, 1, 1.5, and 2.
# This should be an array with the snapshot numbers corresponding with these redshifts.
CMD_SNAPSHOTS=(90 73 61 52 44)

# When the -f option is passed, this script will delete all the existing data products.
# If this is "yes", then doing so will also delete the gadget snapshots.
# IF THE GADGET SNAPSHOTS RETURNED BY get-gadget-snapshot ARE THE ORIGINAL DATA, MAKE SURE THIS IS UNSET.
# OTHERWISE IT WILL DELETE THE ORIGINAL SNAPSHOTS.
CLEAN_GADGET_SNAPS=no
# Should convert-snapshots be run before the rest of the pipeline does?
CONVERT_SNAPSHOTS=yes
# To get CAMELS parameters in the header for the DisPerSE catalogs, this code uses the cosmoastroseed parameter files present in the camels directories.
# If this file is not present, or you simply don't care about having the parameter values in the header, then you can turn this off.
WITH_COSMOASTROSEED=yes
