#!/bin/bash

# Copyright © 2026 Sagan Sutherland
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the “Software”),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the Software
# is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# CAMELS post-processing pipeline
# Sagan Sutherland
#
# There's a lot going on here so I'll do my best to explain it.
# Hopefully that way you can copy this and alter it to fit your own purposes.
# This script will perform a lot of post-processing on simulation outputs.
# Currently the post-processing consists of:
# - Subfind galaxy/halo catalogs
# - Sublink merger trees based on the Subfind catalogs
# - Rockstar galaxy catalogs
# - Consistent Trees merger trees based on the Rockstar catalogs
#
# I've tried to make this script as resumable as possible.
# It'll skip steps that it detects it's already done.
#
# Currently this all assumes you're working on Rusty at the Flatiron Institute.
# If you're working on another cluster, you may have to change which modules get loaded depending on availability, and you'll have to change the paths around.
# You may have to download your own versions of the software being used here.
# This script uses:
# - Arepo (Subfind): https://gitlab.mpcdf.mpg.de/vrs/arepo
# - Sublink: https://ui.adsabs.harvard.edu/abs/2015MNRAS.449...49R/
#            I don't know if this code is publicly available anywhere, but this paper describes the algorithm.
# - Rockstar: https://github.com/jwise77/rockstar-galaxies
# - Consistent Trees: https://bitbucket.org/pbehroozi/consistent-trees/src/main/
# - DisperSE: https://github.com/thierry-sousbie/DisPerSE

PRESET_OPTIONS

##########################
# Shouldn't have to change much beyond this point.
##########################

set -e

function main() {
	if [ "$1" = "check" ]; then
		check
		# Thechnically not necessary since check calls exit
		# But keeping this here conveys the semantics of what happens here locally
		return
	fi

	OPTS=$(getopt -a -n post_processing.sh --options "p:,f" --longoptions "parallel:,force" -- "$@")
	eval set -- "$OPTS"

	while :; do
		case $1 in
		-p | --parallel)
			cpus=$2
			shift 2
			;;
		-f | --force)
			force=
			shift
			;;
		--)
			shift
			break
			;;
		*)
			echo "Unknown option: $1"
			shift
			break
			;;
		esac
	done

	IFS=$'\n'
	mapfile -t ALL_SNAPS < <(sort -g <<<"${ALL_SNAPS[*]}")
	unset IFS
	pwd=$(pwd)
	pwd_without_leading=${pwd#/}
	export SIM_ROOT=${pwd_without_leading////-}
	export HDF5_PLUNGIN_PATH=$(camels-utils hdf5-path)

	# If -f or --force, remove all generated post-processing so we can start fresh
	# If you only want to re-run one particular part of the pipeline, you'll have to manually delete it's corresponding outputs.
	if [ -n "${force+x}" ]; then
		echo "-f passed. Removing all old data products and starting over."
		echo "Press Ctrl-C to cancel now, or forever hold your peace."
		for i in {5..1}; do
			echo -n "$i... "
			sleep 1
		done
		echo 0
		if [ "${CLEAN_GADGET_SNAPS}" = "yes" ] && [ "${CONVERT_SNAPSHOTS}" = "yes" ]; then
			for n in "${ALL_SNAPS[@]}"; do
				rm -f -- "$(get-gadget-snapshot "$n")" &
			done
			rm -f -- "$(get-gadget-ic)" &
		fi
		rm -rf -- "$SUBFIND_OUTPUT" &
		rm -rf -- "$SUBLINK_OUTPUT" &
		rm -rf -- "$SUBLINK_GAL_OUTPUT" &
		rm -rf -- "$ROCKSTAR_OUTPUT" &
		rm -rf -- "$CMD_OUTPUT" &
		rm -rf -- "$PK_OUTPUT" &
		rm -rf -- "$DISPERSE_OUTPUT" &
		wait
	fi

	# If the user does not specify how much parallelism to use, try to determine it from the SLURM variables
	# If the slurm variables do not exist (or are null), defaults to however many cores are available.
	# Note that nproc respects SLURM restrictions on cpus.
	if [ -z "${cpus+x}" ]; then
		cpus=$(nproc)
	fi

	if [ "${CONVERT_SNAPSHOTS}" = "yes" ]; then
		convert-snaps
	fi

	for i in "${ALL_SNAPS[@]}"; do
		if [ ! -e "$(get-gadget-snapshot "$i")" ]; then
			echo "Snapshot $(get-gadget-snapshot "$i") does not exist!"
			missing=
		fi
	done
	if [ -n "${missing+x}" ]; then
		echo "Not proceding without all snapshots present"
		return
	fi

	run-subfind
	run-sublink
	run-rockstar
	run-consistent_tree
	run-cmd
	run-Pk
	run-disperse

	cleanup
}

function run-subfind() {
	read -r Om Ol Ob h box_size < <(camels-utils cosmo-calc $(get-gadget-snapshot "${ALL_SNAPS[0]}"))

	mkdir -p -- "$SUBFIND_OUTPUT"
	cat >"${SUBFIND_OUTPUT}/arepo_subfind_param.txt" <<-EOF
		InitCondFile                                      /dev/null
		Omega0                                            ${Om}
		OmegaLambda                                       ${Ol}
		OmegaBaryon                                       ${Ob}
		HubbleParam                                       ${h}
		BoxSize                                           ${box_size}
		WindEnergyIn1e51erg                               10.3675
		VariableWindVelFactor                             6.76707
		RadioFeedbackFactor                               0.57995
		RadioFeedbackReiorientationFactor                 13.9572
		OutputDir                                         .
		SnapshotFileBase                                  snap
		OutputListFilename                                /dev/null
		ICFormat                                          3
		SnapFormat                                        3
		TimeLimitCPU                                      398000
		CpuTimeBetRestartFile                             3600
		ResubmitOn                                        0
		ResubmitCommand                                   /dev/null
		MaxMemSize                                        6400
		TimeBegin                                         0.0078125
		TimeMax                                           1
		ComovingIntegrationOn                             1
		PeriodicBoundariesOn                              1
		CoolingOn                                         1
		StarformationOn                                   1
		OutputListOn                                      1
		TimeBetSnapshot                                   0
		TimeOfFirstSnapshot                               0
		TimeBetStatistics                                 0.01
		NumFilesPerSnapshot                               1
		NumFilesWrittenInParallel                         16
		TypeOfTimestepCriterion                           0
		ErrTolIntAccuracy                                 0.012
		CourantFac                                        0.3
		MaxSizeTimestep                                   0.005
		MinSizeTimestep                                   2e-09
		CritPhysDensity                                   0
		MaxSfrTimescale                                   2.27
		CritOverDensity                                   57.7
		TempSupernova                                     5.73e+07
		TempClouds                                        1000
		FactorEVP                                         573
		TemperatureThresh                                 0
		InitGasTemp                                       170
		MinGasTemp                                        5
		MinimumDensityOnStartUp                           1e-20
		LimitUBelowThisDensity                            0
		LimitUBelowCertainDensityToThisValue              0
		MinEgySpec                                        0
		TypeOfOpeningCriterion                            1
		ErrTolTheta                                       0.7
		ErrTolForceAcc                                    0.0025
		MultipleDomains                                   8
		TopNodeFactor                                     2.5
		ActivePartFracForNewDomainDecomp                  0.01
		DesNumNgb                                         64
		MaxNumNgbDeviation                                4
		UnitLength_in_cm                                  3.08568e+21
		UnitMass_in_g                                     1.989e+43
		UnitVelocity_in_cm_per_s                          100000
		GravityConstantInternal                           0
		SofteningComovingType0                            $(get-softening 0 comoving)
		SofteningComovingType1                            $(get-softening 1 comoving)
		SofteningComovingType2                            $(get-softening 2 comoving)
		SofteningComovingType3                            $(get-softening 3 comoving)
		SofteningComovingType4                            $(get-softening 4 comoving)
		SofteningMaxPhysType0                             $(get-softening 0 physical)
		SofteningMaxPhysType1                             $(get-softening 1 physical)
		SofteningMaxPhysType2                             $(get-softening 2 physical)
		SofteningMaxPhysType3                             $(get-softening 3 physical)
		SofteningMaxPhysType4                             $(get-softening 4 physical)
		SofteningTypeOfPartType0                          0
		SofteningTypeOfPartType1                          1
		SofteningTypeOfPartType2                          1
		SofteningTypeOfPartType3                          1
		SofteningTypeOfPartType4                          1
		SofteningTypeOfPartType5                          2
		GasSoftFactor                                     2.5
		MinimumComovingHydroSoftening                     0.25
		AdaptiveHydroSofteningSpacing                     1.2
		CellShapingSpeed                                  0.5
		CellMaxAngleFactor                                2.25
		ReferenceGasPartMass                              0
		TargetGasMassFactor                               1
		RefinementCriterion                               1
		DerefinementCriterion                             1
		MHDSeedDir                                        4
		MHDSeedValue                                      1e-14
		TreecoolFile                                      ${AREPO_ROOT}/data/TREECOOL_fg_dec11
		ErrTolThetaSubfind                                0.7
		DesLinkNgb                                        20
		IMF_MinMass_Msun                                  0.1
		IMF_MaxMass_Msun                                  100
		AGB_MassTransferOn                                1
		SNIa_MassTransferOn                               1
		SNII_MassTransferOn                               1
		SNII_MinMass_Msun                                 8
		SNII_MaxMass_Msun                                 100
		SNIa_Rate_TAU                                     0.04
		SNIa_Rate_Norm                                    0.0013
		YieldTablePath                                    ${AREPO_ROOT}/L75TNG_Arepo_GFM_Tables/Yields/
		DesNumNgbEnrichment                               64
		MaxNumNgbDeviationEnrichment                      4
		NSNS_MassTransferOn                               1
		NSNS_MassPerEvent                                 0.05
		NSNS_Rate_TAU                                     0.1
		NSNS_per_SNIa                                     0.001
		ThermalWindFraction                               0.1
		VariableWindSpecMomentum                          0
		WindFreeTravelMaxTimeFactor                       0.025
		WindFreeTravelDensFac                             0.05
		TimeBetOnTheFlyFoF                                1.03
		MinWindVel                                        350
		WindEnergyReductionFactor                         0.25
		WindEnergyReductionMetallicity                    0.002
		WindEnergyReductionExponent                       2
		BlackHoleAccretionFactor                          1
		BlackHoleFeedbackFactor                           0.1
		BlackHoleEddingtonFactor                          1
		SeedBlackHoleMass                                 8e-05
		MinFoFMassForNewSeed                              5
		DesNumNgbBlackHole                                128
		BlackHoleMaxAccretionRadius                       1e+20
		BlackHoleRadiativeEfficiency                      0.2
		QuasarThreshold                                   0.002
		RadioFeedbackMinDensityFactor                     0.01
		BlackHoleCenteringMassMultiplier                  1000
		CoolingTablePath                                  ${AREPO_ROOT}/L75TNG_Arepo_GFM_Tables/Cooling/cooling_metal_AGN_Compton_self_shielding_Rahmati12.hdf5
		MinMetalTemp                                      10000
		PhotometricsTablePath                             ${AREPO_ROOT}/L75TNG_Arepo_GFM_Tables/Photometrics/
		TreecoolFileAGN                                   ${AREPO_ROOT}/data/TREECOOL_AGN
		SelfShieldingDensity                              0.1295
		ObscurationFactor                                 0.3
		ObscurationSlope                                  0.07
		FactorForSofterEQS                                0.3
		TempForSofterEQS                                  10000
		WindDumpFactor                                    0.6
		SelfShieldingFile                                 ${AREPO_ROOT}/data/SelfShielding_Rahmati12
	EOF

	for i in "${ALL_SNAPS[@]}"; do
		# I've seen situations where subfind can spuriously fail, so we end up missing one random catalog.
		# I check each catalog individually rather than checking for the final one so that if one is randomly missing, it can still be recovered without restarting everything.
		if [ -e "$(printf "${SUBFIND_OUTPUT}/groups_%03d.hdf5" "$i")" ]; then
			continue
		fi
		ln -sf "$(realpath "$(get-gadget-snapshot "$i")" --relative-to="${SUBFIND_OUTPUT}")" "$(printf "${SUBFIND_OUTPUT}/snap_%03d.hdf5" "$i")"
		# Make sure Arepo is compiled to output the group-ordered snapshots.
		# Sublink expects these to exist.
		pushd "${SUBFIND_OUTPUT}"
		unset SLURM_JOBID
		mpirun --bind-to none -N "$cpus" "${AREPO_ROOT}/Arepo" ./arepo_subfind_param.txt 3 "$i" 2>&1 | tee "subfind_${i}.log"
		popd
		mv "$(printf "${SUBFIND_OUTPUT}/fof_subhalo_tab_%03d.hdf5" "$i")" "$(printf "${SUBFIND_OUTPUT}/groups_%03d.hdf5" "$i")"
	done
}

function run-sublink() {
	# Subhalos version
	if ! ([ -e "${SUBLINK_OUTPUT}/tree.hdf5" ] &&
		[ -e "${SUBLINK_OUTPUT}/tree_extended.hdf5" ] &&
		[ -e "${SUBLINK_OUTPUT}/offsets/" ]); then

		mkdir -p -- "${SUBLINK_OUTPUT}"

		mkdir -p -- "${SUBLINK_OUTPUT}/subfind"
		n=0
		for i in "${ALL_SNAPS[@]}"; do
			ln -sfr "$(printf "${SUBFIND_OUTPUT}/groups_%03d.hdf5" "$i")" "$(printf "${SUBLINK_OUTPUT}/subfind/fof_subhalo_tab_%03d.hdf5" $n)"
			ln -sfr "$(printf "${SUBFIND_OUTPUT}/snap-groupordered-storeids_%03d.hdf5" "$i")" "$(printf "${SUBLINK_OUTPUT}/subfind/snap_%03d.hdf5" $n)"
			: $((n++))
		done
		((n--))

		# Find descendants
		"${SUBLINK_ROOT}/Descendants/find_descendants" "${SUBLINK_OUTPUT}/subfind" "${SUBLINK_OUTPUT}/sublink_first" 0 $n 0 $n Subhalos first /dev/null
		"${SUBLINK_ROOT}/Descendants/find_descendants" "${SUBLINK_OUTPUT}/subfind" "${SUBLINK_OUTPUT}/sublink_second" 0 $n 0 $n Subhalos second /dev/null

		# Compare descendants
		"${SUBLINK_ROOT}/Descendants/compare_descendants" "${SUBLINK_OUTPUT}/sublink" 0 $n /dev/null

		# Build trees
		"${SUBLINK_ROOT}/SubhaloTrees/build_trees" "${SUBLINK_OUTPUT}/sublink" "${SUBLINK_OUTPUT}/tree" 0 $n /dev/null

		# concatenate, create 'basic' and 'extended' trees, calculate offsets
		# Depending on your system's python environment, you may need to activate your own python to run this.
		python "${SUBLINK_ROOT}/Python/create_columns.py" "${SUBLINK_OUTPUT}/subfind" "${SUBLINK_OUTPUT}" 0 $n 1 0
		python "${SUBLINK_ROOT}/Python/create_extended_trees.py" "${SUBLINK_OUTPUT}/subfind" "${SUBLINK_OUTPUT}" $n
		python "${SUBLINK_ROOT}/Python/concatenate_trees.py" "${SUBLINK_OUTPUT}/tree"
		python "${SUBLINK_ROOT}/Python/concatenate_trees.py" "${SUBLINK_OUTPUT}/tree_extended"
		python "${SUBLINK_ROOT}/Python/compute_offsets.py" "${SUBLINK_OUTPUT}/subfind" "${SUBLINK_OUTPUT}" 0 $n

		# Cleanup
		rm -- "${SUBLINK_OUTPUT}"/sublink_*.hdf5
		rm -r -- "${SUBLINK_OUTPUT}/columns/"
		rm -- "${SUBLINK_OUTPUT}"/tree.*.hdf5
		rm -- "${SUBLINK_OUTPUT}"/tree_extended.*.hdf5
		rm -r -- "${SUBLINK_OUTPUT}/subfind"
	fi

	# Galaxies version
	if ! ([ -e "${SUBLINK_GAL_OUTPUT}/tree.hdf5" ] &&
		[ -e "${SUBLINK_GAL_OUTPUT}/tree_extended.hdf5" ] &&
		[ -e "${SUBLINK_GAL_OUTPUT}/offsets/" ]) &&
		h5ls "$(get-gadget-snapshot "${ALL_SNAPS[0]}")/PartType0" &>/dev/null; then

		mkdir -p -- "${SUBLINK_GAL_OUTPUT}"

		mkdir -p -- "${SUBLINK_GAL_OUTPUT}/subfind"
		n=0
		for i in "${ALL_SNAPS[@]}"; do
			ln -sfr "$(printf "${SUBFIND_OUTPUT}/groups_%03d.hdf5" "$i")" "$(printf "${SUBLINK_GAL_OUTPUT}/subfind/fof_subhalo_tab_%03d.hdf5" $n)"
			ln -sfr "$(printf "${SUBFIND_OUTPUT}/snap-groupordered-storeids_%03d.hdf5" "$i")" "$(printf "${SUBLINK_GAL_OUTPUT}/subfind/snap_%03d.hdf5" $n)"
			: $((n++))
		done
		((n--))

		# Find descendants
		"${SUBLINK_ROOT}/Descendants/find_descendants" "${SUBLINK_GAL_OUTPUT}/subfind" "${SUBLINK_GAL_OUTPUT}/sublink_first" 0 $n 0 $n Galaxies first /dev/null
		"${SUBLINK_ROOT}/Descendants/find_descendants" "${SUBLINK_GAL_OUTPUT}/subfind" "${SUBLINK_GAL_OUTPUT}/sublink_second" 0 $n 0 $n Galaxies second /dev/null

		# Compare descendants
		"${SUBLINK_ROOT}/Descendants/compare_descendants" "${SUBLINK_GAL_OUTPUT}/sublink" 0 $n /dev/null

		# Build trees
		"${SUBLINK_ROOT}/SubhaloTrees/build_trees" "${SUBLINK_GAL_OUTPUT}/sublink" "${SUBLINK_GAL_OUTPUT}/tree" 0 $n /dev/null

		# concatenate, create 'basic' and 'extended' trees, calculate offsets
		# Depending on your system's python environment, you may need to activate your own python to run this.
		python "${SUBLINK_ROOT}/Python/create_columns.py" "${SUBLINK_GAL_OUTPUT}/subfind" "${SUBLINK_GAL_OUTPUT}" 0 $n 1 0
		python "${SUBLINK_ROOT}/Python/create_extended_trees.py" "${SUBLINK_GAL_OUTPUT}/subfind" "${SUBLINK_GAL_OUTPUT}" $n
		python "${SUBLINK_ROOT}/Python/concatenate_trees.py" "${SUBLINK_GAL_OUTPUT}/tree"
		python "${SUBLINK_ROOT}/Python/concatenate_trees.py" "${SUBLINK_GAL_OUTPUT}/tree_extended"
		python "${SUBLINK_ROOT}/Python/compute_offsets.py" "${SUBLINK_GAL_OUTPUT}/subfind" "${SUBLINK_GAL_OUTPUT}" 0 $n

		# Cleanup
		rm -- "${SUBLINK_GAL_OUTPUT}"/sublink_*.hdf5
		rm -r -- "${SUBLINK_GAL_OUTPUT}/columns/"
		rm -- "${SUBLINK_GAL_OUTPUT}"/tree.*.hdf5
		rm -- "${SUBLINK_GAL_OUTPUT}"/tree_extended.*.hdf5
		rm -r -- "${SUBLINK_GAL_OUTPUT}/subfind"
	fi
	rm -f -- "${SUBFIND_OUTPUT}"/snap-groupordered-storeids_*.hdf5
}

function run-rockstar() {
	if [ -e ./rockstar/out_90.list ] || [ -e "${ROCKSTAR_OUTPUT}/trees/tree_0_0_0.dat" ]; then
		return
	fi

	mkdir -p -- "${ROCKSTAR_OUTPUT}"
	# This could probably be done a bit cleaner with a quick python line, but whatever.
	box_size=$(($(h5dump --attribute=/Header/BoxSize --noindex "$(get-gadget-snapshot "${ALL_SNAPS[0]}")" | sed -n "6p") / 1000))

	printf "%03d\n" "${ALL_SNAPS[@]}" >"${ROCKSTAR_OUTPUT}/snaps"

	# TODO: How much does setting the FORCE_RES to the softening length change things?
	# This doesn't matter for most sets, but since SWIMBA SB28 varies the softenings, we should be robust to that variation.
	cat >"${ROCKSTAR_OUTPUT}/rockstar_params.cfg" <<-EOF
		FILE_FORMAT = "AREPO"
		PARTICLE_MASS = 0 # must specify (in Msun/h) for ART or ASCII

		AREPO_LENGTH_CONVERSION = 1e-3
		AREPO_MASS_CONVERSION = 1e+10

		MASS_DEFINITION = "vir" 

		PARALLEL_IO = 1
		PERIODIC = 1

		FORCE_RES = 0.0015

		MIN_HALO_OUTPUT_SIZE = 100

		BOX_SIZE = ${box_size} # Mpc / h

		INBASE = "$(dirname "$(get-gadget-snapshot "${ALL_SNAPS[-1]}")")"
		        FILENAME = "$(basename "$(get-gadget-snapshot "${ALL_SNAPS[-1]}" | sed 's/\(.*\)_.*/\1/')")_<snap>.hdf5"

		SNAPSHOT_NAMES = "${ROCKSTAR_OUTPUT}/snaps

		NUM_BLOCKS = 1

		OUTBASE = "$(realpath "${ROCKSTAR_OUTPUT}")"

		NUM_WRITERS = ${cpus}
		FORK_READERS_FROM_WRITERS = 1
		FORK_PROCESSORS_PER_MACHINE = ${cpus}
	EOF

	# Restart Rockstar if possible
	if [ -e "${ROCKSTAR_OUTPUT}/restart.cfg" ]; then
		params="${ROCKSTAR_OUTPUT}/restart.cfg"
		rm -f -- "${ROCKSTAR_OUTPUT}/auto-rockstar.cfg"
	else
		params="${ROCKSTAR_OUTPUT}/rockstar_params.cfg"
	fi
	# Rockstar works by having one server process and many client processes communicating.
	# Start main process
	"${ROCKSTAR_ROOT}/rockstar-galaxies" -c ${params} &
	while [ ! -f "${ROCKSTAR_OUTPUT}/auto-rockstar.cfg" ]; do
		# The creationg of auto-rockstar.cfg signals that the main process is ready.
		sleep 1
	done
	# Start worker process.
	# This will fork into many workers as necessary.
	"${ROCKSTAR_ROOT}/rockstar-galaxies" -c "${ROCKSTAR_OUTPUT}/auto-rockstar.cfg"

	# Wait for background rockstar process to finish.
	wait

	# Cleanup saved until after consistent trees has been run since it uses many of the rockstar outputs.
}

function run-consistent_tree() {
	if [ -e "${ROCKSTAR_OUTPUT}/trees/tree_0_0_0.dat" ]; then
		return
	fi

	output_dir=$(realpath "${ROCKSTAR_OUTPUT}")

	perl "${ROCKSTAR_ROOT}/scripts/gen_merger_cfg.pl" "${output_dir}/rockstar.cfg"
	# Consistent trees needs to be run inside of its directory, rather than in the directory we exist in.
	# So we cd to there to run these things.
	# This is also why all the paths here and in running rockstar are resolved to absolute paths.
	# Otherwise it'll try to take the relative path from the consistent trees directory, which is wrong.
	pushd "${CONSISTENT_TREE_ROOT}"
	while
		perl "${CONSISTENT_TREE_ROOT}/do_merger_tree.pl" "${output_dir}/outputs/merger_tree.cfg" 2>&1 | tee "${output_dir}/mergers.log"
		grep "Error: too few halos at scale factor" "${output_dir}/mergers.log"
	do
		# Remove the first scale factor and try again.
		tail -n+2 -- "${output_dir}/outputs/scales.txt" >"${output_dir}/outputs/scales2.txt"
		mv -- "${output_dir}/outputs/scales2.txt" "${output_dir}/outputs/scales.txt"
	done
	perl "${CONSISTENT_TREE_ROOT}/halo_trees_to_catalog.pl" "${output_dir}/outputs/merger_tree.cfg"
	popd
	if [ ! -e "${output_dir}/trees/tree_0_0_0.dat" ]; then
		echo "no ${output_dir}/trees/tree_0_0_0.dat, probably failed"
		return
	else
		echo "${output_dir}/trees/tree_0_0_0.dat exists, probably succeeded"
	fi

	# Cleanup once everything's done.
	rm -f -- "${output_dir}"/halos*.{ascii,bin}
	rm -f -- "${output_dir}"/out_*.list
	rm -f -- "${output_dir}"/{auto-rockstar,restart,rockstar}.cfg
	rm -rf -- "${output_dir}/outputs/" "${output_dir}/profiling/"
}

function run-cmd() {
	# CAMELS Multifield Dataset
	# make-CMD is included in my postprocessing package.
	# It's mainly based on Paco's CMD code, but I've remade it from scratch to hopefully make it a bit easier to extend.
	# I've also made it parallelized without using mpi, since I've historically had consistency issues with mpi4py.

	snapshot="$(get-gadget-snapshot "${ALL_SNAPS[0]}")"
	if h5ls "${snapshot}/PartType0" &>/dev/null; then
		target_2d=12
		target_3d=180
	else
		target_2d=1
		target_3d=15
	fi

	mkdir -p -- "${CMD_OUTPUT}/2D_maps"
	mkdir -p -- "${CMD_OUTPUT}/3D_grids"
	if [ "$(ls -- "${CMD_OUTPUT}"/2D_maps/ | wc -l)" -lt "$target_2d" ]; then
		OMP_NUM_THREADS="$cpus" make-CMD "$(get-gadget-snapshot "${ALL_SNAPS[-1]}")" --parallel "$cpus" --target "${CMD_OUTPUT}/2D_maps" --grid 256 --2d
	fi
	if [ "$(ls -- "${CMD_OUTPUT}"/3D_grids/ | wc -l)" -lt "$target_3d" ]; then
		IFS=$'\n' # Split the snapshot names into different arguments
		OMP_NUM_THREADS="$cpus" make-CMD $(for n in "${CMD_SNAPSHOTS[@]}"; do
			get-gadget-snapshot "${n}"
			echo
		done) --parallel "$cpus" --target "${CMD_OUTPUT}/3D_grids" --grid 128 --3d
		OMP_NUM_THREADS="$cpus" make-CMD $(for n in "${CMD_SNAPSHOTS[@]}"; do
			get-gadget-snapshot "${n}"
			echo
		done) --parallel "$cpus" --target "${CMD_OUTPUT}/3D_grids" --grid 256 --3d
		OMP_NUM_THREADS="$cpus" make-CMD $(for n in "${CMD_SNAPSHOTS[@]}"; do
			get-gadget-snapshot "${n}"
			echo
		done) --parallel "$cpus" --target "${CMD_OUTPUT}/3D_grids" --grid 512 --3d
		unset IFS
	fi
}

function run-Pk() {
	mkdir -p -- "${PK_OUTPUT}"
	pk_files=("${PK_OUTPUT}"/*)
	if [ ${#pk_files[@]} -ge 460 ]; then
		return
	fi
	make-Pk "$(get-gadget-snapshot "${ALL_SNAPS[-1]}" | sed 's/\(.*\)_.*/\1/')"_*.hdf5 --parallel "${cpus}" --target "${PK_OUTPUT}"
	ic="$(get-gadget-ic)"
	if [ -e "$ic" ]; then
		make-Pk "$ic" --parallel "$cpus" --target "${PK_OUTPUT}"
	fi
}

function run-disperse {
	grid=256
	cut=500
	sigma=2
	snap=${ALL_SNAPS[-1]}
	if [ "$cut" -ge 1000 ]; then
		cut_str=$(printf "%.0e" "$cut")
	else
		cut_str=$(printf "%d" "$cut")
	fi
	sigma_str=$(printf "%02d" "$sigma")
	snap_str=$(printf "%03d" "$snap")

	if [ -e "${DISPERSE_OUTPUT}/massgrid-disperse_G-${grid}_S-${sigma_str}_c${cut_str}-output_file_${snap_str}.hdf5" ]; then
		return
	fi

	mkdir -p -- "${DISPERSE_OUTPUT}"

	read -r Om Ol Ok h w box_size < <(python -c "import h5py; f=h5py.File('$(get-gadget-snapshot "$snap")'); h=f['Header'].attrs; print(h['Omega0'], h['OmegaLambda'], 1 - h['Omega0'] - h['OmegaLambda'], h['HubbleParam'], '-1')")

	# make-mesh is included as part of my postprocessing package
	make-mesh "$(get-gadget-snapshot "$snap")" --parallel "$cpus" --target "${DISPERSE_OUTPUT}" --grid "$grid" --sigma "$sigma"
	"${DISPERSE_ROOT}/fieldconv" "${DISPERSE_OUTPUT}/masscubegrid-G-${grid}_S-${sigma_str}_${snap_str}.fits" -cosmo "$Om" "$Ol" "$Ok" "$h" "$w" -to NDfield -info -outName "masscubegrid-G-${grid}_S-${sigma_str}_${snap_str}" -outDir "${DISPERSE_OUTPUT}"
	if [ -e "${DISPERSE_OUTPUT}/masscubegrid-G-256_S-${sigma_str}_${snap_str}.MSC" ]; then
		"${DISPERSE_ROOT}/mse" "${DISPERSE_OUTPUT}/masscubegrid-G-${grid}_S-${sigma_str}_${snap_str}.ND" -cut $cut -upSkl -manifolds -forceLoops -periodicity 111 -nthreads "$cpus" -outName "masscubegrid-G-${grid}_S-${sigma_str}_${snap_str}" -outDir "${DISPERSE_OUTPUT}" -loadMSC "${DISPERSE_OUTPUT}/masscubegrid-G-256_S-02_${snap_str}.MSC"
	else
		"${DISPERSE_ROOT}/mse" "${DISPERSE_OUTPUT}/masscubegrid-G-${grid}_S-${sigma_str}_${snap_str}.ND" -cut $cut -upSkl -manifolds -forceLoops -periodicity 111 -nthreads "$cpus" -outName "masscubegrid-G-${grid}_S-${sigma_str}_${snap_str}" -outDir "${DISPERSE_OUTPUT}"
	fi
	"${DISPERSE_ROOT}/skelconv" "${DISPERSE_OUTPUT}/masscubegrid-G-${grid}_S-${sigma_str}_${snap_str}_c${cut_str}.up.NDskl" -cosmo "$Om" "$Ol" "$Ok" "$h" "$w" -breakdown -smooth 1 -to segs_ascii -outName "masscubegrid-G-${grid}_S-${sigma_str}_${snap_str}_c${cut_str}" -outDir "${DISPERSE_OUTPUT}"
	"${DISPERSE_ROOT}/skelconv" "${DISPERSE_OUTPUT}/masscubegrid-G-${grid}_S-${sigma_str}_${snap_str}_c${cut_str}.up.NDskl" -cosmo "$Om" "$Ol" "$Ok" "$h" "$w" -breakdown -smooth 1 -to crits_ascii -outName "masscubegrid-G-${grid}_S-${sigma_str}_${snap_str}_c${cut_str}" -outDir "${DISPERSE_OUTPUT}"
	"${DISPERSE_ROOT}/skelconv" "${DISPERSE_OUTPUT}/masscubegrid-G-${grid}_S-${sigma_str}_${snap_str}_c${cut_str}.up.NDskl" -cosmo "$Om" "$Ol" "$Ok" "$h" "$w" -breakdown -smooth 1 -to NDskl_ascii -outName "masscubegrid-G-${grid}_S-${sigma_str}_${snap_str}_c${cut_str}" -outDir "${DISPERSE_OUTPUT}"
	if [ "${WITH_COSMOASTROSEED}" = "yes" ]; then
		# disperse2hdf5 is part of postprocessing package
		disperse2hdf5 "${DISPERSE_OUTPUT}/masscubegrid-G-${grid}_S-${sigma_str}_${snap_str}_c${cut_str}.BRK.S001.a.NDskl" "${SUBFIND_OUTPUT}/groups_${snap_str}.hdf5" ../CosmoAstroSeed*.txt --grid $grid --target "${DISPERSE_OUTPUT}"
	else
		disperse2hdf5 "${DISPERSE_OUTPUT}/masscubegrid-G-${grid}_S-${sigma_str}_${snap_str}_c${cut_str}.BRK.S001.a.NDskl" "${SUBFIND_OUTPUT}/groups_${snap_str}.hdf5" --grid $grid --target "${DISPERSE_OUTPUT}"
	fi

	# Cleanup
	find "${DISPERSE_OUTPUT}" -mindepth 1 -not -name '*.hdf5' -delete
	mv "${DISPERSE_OUTPUT}/masscubegrid-G-${grid}_S-${sigma_str}_${snap_str}_c${cut_str}-output_file_${snap_str}.hdf5" "${DISPERSE_OUTPUT}/massgrid-disperse_G-${grid}_S-${sigma_str}_c${cut_str}-output_file_${snap_str}.hdf5"
}

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

function check() {
	n_snap=${#ALL_SNAPS[@]}
	snapshot="$(get-gadget-snapshot "${ALL_SNAPS[0]}")"
	if h5ls "${snapshot}/PartType0" &>/dev/null; then
		n_body="n"
	else
		n_body="y"
	fi

	result=0
	ensure_count "subfind/groups_*" $n_snap
	ensure_exists sublink/tree.hdf5 sublink/tree_extended.hdf5
	ensure_count "sublink/offsets/offsets_*.hdf5" $n_snap
	ensure_exists sublink_gal/tree.hdf5 sublink_gal/tree_extended.hdf5
	ensure_count "sublink_gal/offsets/offsets_*.hdf5" $n_snap
	ensure_exists rockstar/trees/tree_0_0_0.dat rockstar/trees/locations.dat rockstar/trees/forests.list
	ensure_count rockstar/hlists/ $(($n_snap - 20)) # This isn't an exact count since it'll skip snapshots with not enough halos
	if [ $n_body = "n" ]; then
		ensure_count CMD/2D_maps 12
		ensure_count CMD/3D_grids $((3 * 5 * 12))
		ensure_count Pk/ $((($n_snap + 1) * 5))
	else
		ensure_count CMD/2D_maps 1
		ensure_count CMD/3D_grids $((3 * 5))
		ensure_count Pk/ $(($n_snap + 1))
	fi
	ensure_exists DisPerSE/massgrid-disperse_G-256_S-02_c500-output_file_090.hdf5

	if [ "$result" -eq 0 ]; then
		exit 0
	else
		exit 1
	fi
}

main "$@"
