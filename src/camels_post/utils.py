from pathlib import Path
import argparse


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(required=True)

    hdf5_path = subparsers.add_parser(
        "hdf5-path",
        help="Print the HDF5 plugin path from hdf5plugin",
    )
    hdf5_path.set_defaults(func=print_hdf5)

    cosmo = subparsers.add_parser(
        "cosmo-calc",
        help="Determine the cosmology of a simulation",
    )
    cosmo.add_argument("snapshot")
    cosmo.set_defaults(func=calculate_cosmology)

    args = parser.parse_args()
    args.func(args)


def print_hdf5(_):
    import hdf5plugin

    print(hdf5plugin.PLUGIN_PATH)


def calculate_cosmology(args):
    import h5py
    import hdf5plugin as _  # noqa
    import numpy as np

    with h5py.File(args.snapshot) as snap:
        h = snap["Header"].attrs
        M = 0.0
        Mb = 0.0
        for p in (1, 2):
            pt = "PartType" + str(p)
            if pt not in snap:
                continue
            if "Masses" in snap[pt]:
                M += np.sum(snap[pt]["Masses"][:])
            else:
                M += h["MassTable"][p] * h["NumPart_Total"][p]
        for p in (0, 3, 4, 5, 6):
            pt = "PartType" + str(p)
            if pt not in snap:
                continue
            if "Masses" in snap[pt]:
                Mb += np.sum(snap[pt]["Masses"][:])
            else:
                Mb += h["MassTable"][p] * h["NumPart_Total"][p]
        M += Mb
        print(
            "%6f %6f %6f %6f %d"
            % (
                h["Omega0"],
                h["OmegaLambda"],
                h["Omega0"] * Mb / M,
                h["HubbleParam"],
                round(h["BoxSize"]),
            )
        )
