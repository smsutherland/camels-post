"""
Make-mesh is largely based on code provided by Jonathan Mercedes Feliz.
"""

import argparse
from pathlib import Path

import h5py
import hdf5plugin as _  # noqa
import joblib
import numpy as np
from astropy.io import fits
from meshoid import Meshoid
from scipy import signal


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "snapshots",
        nargs="+",
        help="gadget-hdf5 style snapshot(s) to deposit onto a grid",
        type=Path,
    )
    parser.add_argument("--target", type=Path, default=Path.cwd())
    parser.add_argument("--grid", type=int, default=256)
    parser.add_argument("--sigma", type=int, default=2)
    parser.add_argument("--parallel", "-p", type=int, default=joblib.cpu_count())
    parser.add_argument("--gas", action="store_true")

    args = parser.parse_args()
    snapshots: list[Path] = args.snapshots
    target: Path = args.target
    grid: int = args.grid
    sigma: int = args.sigma
    parallel: int = args.parallel
    gas: bool = args.gas

    ptype = 0 if gas else 1

    target.mkdir(parents=True, exist_ok=True)

    parallel = min(parallel, len(snapshots))
    joblib.parallel.Parallel(n_jobs=parallel)(
        joblib.delayed(make_mesh)(
            snap, grid=grid, sigma=sigma, target=target, ptype=ptype
        )
        for snap in snapshots
    )


def make_mesh(snap: Path, grid: int, sigma: int, target: Path, ptype: int = 1):
    with h5py.File(snap) as f:
        coordinates: np.ndarray = f[f"PartType{ptype}/Coordinates"][:]
        if "Masses" in f[f"PartType{ptype}"]:
            masses: np.ndarray = f[f"PartType{ptype}/Masses"][:] * 1e10
        else:
            masses = np.full(
                coordinates.shape[0], f["Header"].attrs["MassTable"][ptype] * 1e10
            )
        boxsize: float = f["Header"].attrs["BoxSize"]
        coordinates %= boxsize
    mesh = Meshoid(
        coordinates,
        m=masses,
        des_ngb=32,
        boxsize=boxsize,
        verbose=True,
    )
    massgrid = mesh.DepositToGrid(
        masses,
        res=grid,
        size=boxsize,
        center=np.array(
            [
                (boxsize - (boxsize / grid)) / 2,
            ]
            * 3
        ),
    )

    snapnum = snap.stem[-3:]

    if sigma > 0:
        kernel = gaussian_kernel_3D(sigma, 3)

        # TODO: make a periodic convolution?
        filtered_data = signal.convolve(massgrid, kernel, mode="same")

        filtered_hdu = fits.PrimaryHDU(filtered_data)
        filtered_hdul = fits.HDUList([filtered_hdu])
        filtered_hdul.writeto(
            target / f"masscubegrid-G-{grid}_S-{sigma:02d}_{snapnum}.fits",
            overwrite=True,
        )
    else:
        hdu = fits.PrimaryHDU(massgrid)
        hdul = fits.HDUList([hdu])
        hdul.writeto(
            target / f"masscubegrid-G-{grid}_{snapnum:03d}.fits",
            overwrite=True,
        )


def gaussian_kernel_3D(sigma: float, sigma_level: float):
    sigma_extent = sigma_level * sigma
    extent = np.arange(-sigma_extent, sigma_extent + 1, 1)

    xx, yy, zz = np.meshgrid(extent, extent, extent)
    kernel = (1.0 / (sigma**3 * (2 * np.pi) ** (3 / 2))) * np.exp(
        -(xx**2 + yy**2 + zz**2) / (2 * sigma**2)
    )

    return kernel


if __name__ == "__main__":
    main()
