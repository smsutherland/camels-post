camels-post
-----------
`camels-post` is a all-in-one package for post-processing simulations for the
[CAMELS project](https://camels.readthedocs.io). The package is designed to work
out of the box, though it has currently only been tested on the Flatiron
Institute's rusty and popeye clusters. The package supports both hydrodynamic
and N-body simulations. If you encounter any problems using `camels-post`,
please file an issue on github. If you use `camels-post` to generate data for a
publication, I ask that you please note this in an acknowledgments section.

### Quick Navigation
  * [Data Products](#data_products)
  * [Installation](#installation)
  * [Usage](#usage)
  * [Configuration](#configuration)

### Data Products

#### Subfind halo/subhalo catalogs

The subfind algorithm, implemented in the
[Arepo](https://gitlab.mpcdf.mpg.de/vrs/arepo) simulation code, is used to
identify substructure within snapshots. Outputs are named
`groups_<snap_num>.hdf5` and are put in the directory specified by
[SUBFIND_OUTPUT](#_output). Each hdf5 output contains both friends-of-friends
halos and subfind subhalos. A good description of the data can be found on the
[IllustrisTNG
website](https://www.tng-project.org/data/docs/specifications/#sec2). Arepo is
not shipped with this package. You must provide your own copy, and point to it's
home directory via the [AREPO_ROOT](#_root) config parameter. Based on my own
testing, the public version of Arepo can kinda work, but doesn't produce all the
expected outputs in its files. If you have access to it, I recommend using the
build of Arepo used by the CAMELS project for the IllustrisTNG suite.

#### Sublink merger trees

#### Rockstar + Consistent Trees merger trees

#### DisPerSE cosmic web filaments

#### Camels Multifield Dataset (CMD) maps and grids

#### Power spectra

### Installation
`camels-post` depends on [Eigen](https://libeigen.gitlab.io/) (Through the
[voxelize](https://github.com/leanderthiele/voxelize) library).
It also requires a python installation with development headers Cython can use.
`camels-post` also requires an MPI installation in order to run some parts of
the pipeline. `camels-post` is designed for use with
[disBatch](https://github.com/flatironinstitute/disBatch), but it is not
strictly required.

Installation can be achieved via pip:

```
pip install 'git+https://github.com/smsutherland/camels-post'
```

Or via uv:

```
uv tool install 'git+https://github.com/smsutherland/camels-post'
```

Or via pipx:

```
pipx install 'git+https://github.com/smsutherland/camels-post'
```

#### A note on modules

### Usage
The main tool you will interact with is called `camels-post`.

```
usage: camels-post [-h] [--target TARGET] [--cpus-per-task CPUS_PER_TASK]
                   [--ntasks NTASKS] [--force] [--preset {swimba,gadget-nbody}]
                   dirs [dirs ...]

positional arguments:
  dirs                  Where to look for simulations to post-process

options:
  -h, --help            show this help message and exit
  --target TARGET       Directory to put sbatch scripts
  --cpus-per-task CPUS_PER_TASK, -c CPUS_PER_TASK
                        How many cpus should each post-processing task use?
  --ntasks NTASKS, -n NTASKS
                        How many post-processing tasks should run at once?
  --force, -f           Make post-processing jobs start over
  --preset {swimba,gadget-nbody}
                        Which preset post-processing configuration to use?
```

For instance, if you want to process your CV set, run
```
camels-post path/to/your/suite/CV/CV_*
```
The following files will be added to your current directory:
  * post_process.sh
  * tasks
  * check_tasks
  * job.sh

post_process.sh contains the bash script to performing processing on a single
run. I highly recommend looking at this script. Editing this script is how
you configure `camels-post`. See [Configuration](#configuration). tasks
contains a disBatch task list to process each run. check_tasks contains a
disBatch task list to verify the outputs of the pipeline. Note that this step
only verifies the existence of outputs. This makes no guarantee about the
correctness of outputs. job.sh contains a
[sbatch](https://slurm.schedmd.com/overview.html) script for launching the
post-processing tasks.

### Adding a preset configuration

Preset configurations are located [here](src/camels_post/presets/). Currently
the only presets available are for SWIMBA and Gadget-nbody. A preset consists
of all the options and functions which the pipeline relies on to work for your
specific simulation setup (see [Configuration](#configuration)). To add a new
preset, copy an existing one and change the settings as you see fit. You will
need to add the name of your preset as an option in
[postprocess.py](src/camels_post/postprocess.py#L42).

### Configuration

The following functions are defined in a post-processing configuration.

#### convert-snaps()

Should perform any conversions necessary
such that your snapshots *and* initial conditions are in a Gadget hdf5 format.
If no conversions are needed, set [CONVERT_SNAPSHOTS](#convert_snapshots) to
"no".

#### get-gadget-ic()

Should print a path to the initial conditions. ICs and snapshots are expected
to exist as a single hdf5 file with Gadget's unit conventions (10^10 Msun/h,
kpc/h, etc.). `combine-IC` is a utility shipped with this package to take a
Gadget multi-file IC (either format 1 or hdf5) and turn it into a single hdf5
file.

#### get-gadget-snapshot(num)

Print the path to a gadget hdf5 snapshot. Num ($1) determines the snapshot
number to print (see [ALL_SNAPS](#all_snaps)).

#### cleanup()

This function is called after the pipeline is successfully completed. If you
created any temporary files (such as in [convert-snaps](#convert-snaps)), you
should delete them here.

#### get-softening(ptype, physical/comoving)

[Subfind](#subfind-halosubhalo-catalogs) depends on the gravitational softening
parameters matching what the simulation was run with. This function should print
the softening length for the given particle type in kpc/h, either the comoving
length, or the maximum physical length. Parameter 1 is the numerical particle
type (0, 1, 4, 5). Parameter 2 is either "comoving" or "physical".

The following variables are also defined. Remember that this config is run as
part of a larger bash script, so these variables can have values dynamically
determined at runtime!

#### ALL_SNAPS

This parameter should be an array containing all the snap numbers your
simulation produces. The numbers do not need to be zero-padded, but not
zero-padding them will mean you need to handle any zero-padding in
[get-gadget-snapshot](#get-gadget-snapshotnum)

#### *_OUTPUT

The various *_OUTPUT variables determine where to put the outputs for each step
of the process. Currently, the steps to specify outputs for are SUBFIND,
SUBLINK, SUBLINK_GAL, ROCKSTAR, CMD, PK, and DISPERSE.

#### *_ROOT

Some of the data products come from tools provided in this package (CMD, power
spectra). Others depend on external tools. `camels-post` uses the *_ROOT
variables to know where these tools are installed. Currently, the tools to
specify roots for are AREPO (subfind), SUBLINK, ROCKSTAR, CONSISTENT_TREE,
DISPERSE. Note that Sublink has multiple versions (multi vs single file
snapshots, 4 vs 8 byte particle IDs). Other parts of the pipeline require
single file snapshots ([so convert as necessary](#convert-snaps)), but make sure
to use the correct ID length.

#### *_PTYPES

This is a WIP feature currently. BARYON_PTYPES should be an array with the
particle types which are baryons. For a hydrodynamic run this should be 0, 4,
and 5 for compatability with all the tools. For dark matter only runs, this
array should be empty. DM_PTYPES should be an array of just 1. This option will
eventually be removed.

#### CMD_SNAPSHOTS

In the CAMELS Multifield Dataset, we make 3D grids for z = 0, 0.5, 1, 1.5, and
2. This should be an array with the snapshot numbers corresponding with these
redshifts.

#### CLEAN_GADGET_SNAPS

When `--force` is passed to the post-processing script, it will remove all its
existing data products and start over from scratch. If this option is set to
"yes", the snapshots will be removed as well. **Only set this option to "yes" if
your snapshots are converted. If [get-gadget-snapshot](#get-gadget-snapshotnum)
points to your original data products, turn this option off.**

#### CONVERT_SNAPSHOTS

If set to "yes", the pipeline will call [convert-snaps](#convert-snaps) at the
start.

#### WITH_COSMOASTROSEED

DisPerSE outputs include the parameters that go in to a given CAMELS run. The
parameters are given by the CosmoAstroSeed*.txt file in the parent directory. If
no such file exists, turn this option off by setting it to "no".

### To-Dos
- Make the preset list populate automatically based on the contents of the
  preset directory.
- Move swift2gadget to this package
- 50Mpc
- Allow snapshots to not all be in the form path/to/snap_<snapNumber>.hdf5
- Either make all the parts use BARYON_PTYPES and DM_PTYPES, or enforce
  everything being like SIMBA. Currently the pipeline breaks if, for example,
  a snapshot has background DM type 2.
- restart_job.sh
- Verify configuration validity at start of script.
- Automatically determine CMD snapshots from redshifts in the snapshots.
- Document snapshot format minimum requirements.
- Improve my tooling to use multi-file snapshots.
- Combine-IC should make a virtual snapshot if combining hdf5 snapshots.
- Rename combine-IC to combine-file.
- Try to reconstruct groupordered snapshots based on subfind output.
- Option to control which steps run.
- Rename outputs to better match names in CAMELS.
- Check how the subfind differ between our arepo and the public version.
- Make modules more generalisable.
- specify MPI dependency
- Determine minimum Arepo build
