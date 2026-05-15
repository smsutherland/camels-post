module purge
module load modules/2.3
module load gcc
module load openmpi
module load hdf5/mpi
module load fftw/mpi
module load openblas/threaded
module load flexiblas
module load cfitsio
module load gsl
module load gmp
module load eigen
module load python

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
ROCKSTAR_ROOT=/mnt/home/ssutherland/codes/rockstar-galaxies_dm/
CONSISTENT_TREE_ROOT=/mnt/sdceph/users/camels/Codes/consistent-trees/
DISPERSE_ROOT=/mnt/home/ssutherland/codes/DisPerSE/bin/
