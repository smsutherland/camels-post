import argparse
import os
from dataclasses import dataclass
from pathlib import Path

import h5py
import joblib
import MAS_library.MAS_library as MASL
import numpy as np
from scipy.spatial import KDTree
from voxelize import Voxelize


@dataclass(slots=True, kw_only=True)
class SnapshotData:
    gas_position: np.ndarray  # Mpc / h
    gas_velocity: np.ndarray  # km / s
    gas_mass: np.ndarray  # h Msun / Mpc^2
    gas_temperature: np.ndarray  # K
    gas_pressure: np.ndarray  # h^2 Msun (km / s)^2 / kpc^3
    gas_metallicity: np.ndarray  # dimensionless
    gas_neutral_H: np.ndarray  # h Msun / Mpc^2
    gas_electron_density: np.ndarray  # h^3 / cm^3 / Mpc
    gas_radius: np.ndarray  # Mpc / h
    gas_magnesium: np.ndarray  # dimensionless
    gas_iron: np.ndarray  # dimensionless

    dm_position: np.ndarray  # Mpc / h
    dm_mass: np.ndarray  # h Msun / Mpc^2
    dm_velocity: np.ndarray  # km / s
    dm_radius: np.ndarray  # Mpc / h

    star_position: np.ndarray  # Mpc / h
    star_mass: np.ndarray  # h Msun / Mpc^2
    star_radius: np.ndarray  # Mpc / h

    bh_position: np.ndarray  # Mpc / h
    bh_mass: np.ndarray  # Mpc / h
    bh_radius: np.ndarray  # Mpc / h

    box_size: float  # Mpc / h
    redshift: float
    has_parts: np.ndarray


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "snapshots",
        nargs="+",
        help="A gadget-hdf5 format snapshot(s) which will be imaged.",
        type=Path,
    )
    parser.add_argument("--target", type=Path, default=Path.cwd())
    parser.add_argument(
        "-p",
        "--parallel",
        help="How many parallel tasks to run.",
        type=int,
        default=joblib.cpu_count(),
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("--grid", type=int, default=256)
    parser.add_argument("--r-divisions", type=int, default=20)
    parser.add_argument("--tracers", type=int, default=1000)
    parser.add_argument("--splits", type=int, default=5)
    parser.add_argument("--2d", action="store_true", help="Only do 2d maps")
    parser.add_argument("--3d", action="store_true", help="Only do 3d maps")

    args = parser.parse_args()
    snaps: list[Path] = args.snapshots
    target_dir: Path = args.target
    parallelism: int = args.parallel
    verbose: bool = args.verbose
    grid: int = args.grid
    r_divisions: int = args.r_divisions
    tracers: int = args.tracers
    splits: int = args.splits
    run_2d: bool = getattr(args, "2d")
    run_3d: bool = getattr(args, "3d")

    target_dir.mkdir(parents=True, exist_ok=True)
    snap_data = [load_snap(snap, parallelism) for snap in snaps]

    if run_2d or not run_3d:
        n_jobs = min(parallelism, 3 * splits * len(snaps))
        results = joblib.parallel.Parallel(
            n_jobs=n_jobs,
            verbose=int(verbose) * 10,
            return_as="generator",
        )(
            joblib.parallel.delayed(make_images)(
                axis,
                snap.box_size * slice / splits,
                snap.box_size * (slice + 1) / splits,
                snap,
                verbose=verbose,
                grid=grid,
                r_divisions=r_divisions,
                tracers=tracers,
            )
            for snap in snap_data
            for axis in range(3)
            for slice in range(splits)
        )

        fields = [
            "Mgas",
            "Vgas",
            "Mcdm",
            "Vcdm",
            "Mstar",
            "Mtot",
            "T",
            "Z",
            "P",
            "HI",
            "ne",
            "MgFe",
        ]

        for _ in range(len(snaps)):
            full_results = {
                field: np.empty((3 * splits, grid, grid), dtype=np.float64)
                for field in fields
            }
            z = 0
            for i in range(3 * splits):
                result = next(results)
                for f in fields:
                    if f in result:
                        full_results[f][i] = result[f]
                    else:
                        if f in full_results:
                            full_results.pop(f)
                z = result["z"]

            # This does not include the code name or the set name.
            # Those should be added in post when combining files together.
            suffix = f"z={z:.2f}.npy"
            paths = {field: target_dir / f"Maps_{field}_{suffix}" for field in fields}

            for field, path in paths.items():
                if field in full_results:
                    np.save(path, full_results[field])

    if run_3d or not run_2d:
        os.environ["OMP_NUM_THREADS"] = str(parallelism)
        results = (
            make_grids(
                snap,
                verbose=verbose,
                grid=grid,
            )
            for snap in snap_data
        )

        fields = [
            "Mgas",
            "Vgas",
            "Mcdm",
            "Vcdm",
            "Mstar",
            "Mtot",
            "T",
            "Z",
            "P",
            "HI",
            "ne",
            "MgFe",
        ]

        for snap in snaps:
            result = next(results)
            z = result["z"]

            # This does not include the code name or the set name.
            # Those should be added in post when combining files together.
            suffix = f"{grid}_z={z:.2f}.npy"
            paths = {field: target_dir / f"Grids_{field}_{suffix}" for field in fields}

            for field, path in paths.items():
                if field in result:
                    np.save(path, result[field])

    # shutil.rmtree(target_dir / "mmap")


def make_images(
    axis: int,
    pos_min: float,
    pos_max: float,
    data: SnapshotData,
    *,
    verbose: bool = False,
    grid: int = 256,
    r_divisions: int = 20,
    tracers: int = 1000,
):
    axis_coord = data.gas_position[:, axis]
    indices_gas = np.nonzero((axis_coord >= pos_min) & (axis_coord < pos_max))[0]
    axis_coord = data.dm_position[:, axis]
    indices_dm = np.nonzero((axis_coord >= pos_min) & (axis_coord < pos_max))[0]
    axis_coord = data.star_position[:, axis]
    indices_star = np.nonzero((axis_coord >= pos_min) & (axis_coord < pos_max))[0]
    axis_coord = data.bh_position[:, axis]
    indices_bh = np.nonzero((axis_coord >= pos_min) & (axis_coord < pos_max))[0]

    axis_x = (axis + 1) % 3
    axis_y = (axis + 2) % 3

    masked_data = SnapshotData(
        gas_position=data.gas_position[:, [axis_x, axis_y]][indices_gas],
        gas_velocity=data.gas_velocity[indices_gas],
        gas_mass=data.gas_mass[indices_gas],
        gas_temperature=data.gas_temperature[indices_gas],
        gas_pressure=data.gas_pressure[indices_gas],
        gas_metallicity=data.gas_metallicity[indices_gas],
        gas_neutral_H=data.gas_neutral_H[indices_gas],
        gas_electron_density=data.gas_electron_density[indices_gas],
        gas_radius=data.gas_radius[indices_gas],
        gas_magnesium=data.gas_magnesium[indices_gas],
        gas_iron=data.gas_iron[indices_gas],
        dm_position=data.dm_position[:, [axis_x, axis_y]][indices_dm],
        dm_mass=data.dm_mass[indices_dm],
        dm_velocity=data.dm_velocity[indices_dm],
        dm_radius=data.dm_radius[indices_dm],
        star_position=data.star_position[:, [axis_x, axis_y]][indices_star],
        star_mass=data.star_mass[indices_star],
        star_radius=data.star_radius[indices_star],
        bh_position=data.bh_position[:, [axis_x, axis_y]][indices_bh],
        bh_mass=data.bh_mass[indices_bh],
        bh_radius=data.bh_radius[indices_bh],
        box_size=data.box_size,
        redshift=data.redshift,
        has_parts=data.has_parts,
    )

    box_size = masked_data.box_size

    fields = {}

    if masked_data.has_parts[0]:
        fields |= {
            "T": (
                masked_data.gas_position,
                masked_data.gas_radius,
                masked_data.gas_temperature * masked_data.gas_mass,
            ),
            "Z": (
                masked_data.gas_position,
                masked_data.gas_radius,
                masked_data.gas_metallicity * masked_data.gas_mass,
            ),
            "Vgas": (
                masked_data.gas_position,
                masked_data.gas_radius,
                masked_data.gas_velocity * masked_data.gas_mass,
            ),
            "P": (
                masked_data.gas_position,
                masked_data.gas_radius,
                masked_data.gas_pressure * masked_data.gas_mass,
            ),
            "Mgas": (
                masked_data.gas_position,
                masked_data.gas_radius,
                masked_data.gas_mass,
            ),
            "HI": (
                masked_data.gas_position,
                masked_data.gas_radius,
                masked_data.gas_neutral_H,
            ),
            "ne": (
                masked_data.gas_position,
                masked_data.gas_radius,
                masked_data.gas_electron_density,
            ),
            "Mg": (
                masked_data.gas_position,
                masked_data.gas_radius,
                masked_data.gas_magnesium,
            ),
            "Fe": (
                masked_data.gas_position,
                masked_data.gas_radius,
                masked_data.gas_iron,
            ),
        }
    if masked_data.has_parts[1] and masked_data.has_parts.sum() > 1:
        # Only include this if we're in a hydro run.
        # nbody runs get just m_tot
        fields |= {
            "Mcdm": (
                masked_data.dm_position,
                masked_data.dm_radius,
                masked_data.dm_mass,
            ),
            "Vcdm": (
                masked_data.dm_position,
                masked_data.dm_radius,
                masked_data.dm_velocity * masked_data.dm_mass,
            ),
        }
    if masked_data.has_parts[4]:
        fields |= {
            "Mstar": (
                masked_data.star_position,
                masked_data.star_radius,
                masked_data.star_mass,
            ),
        }

    results = {}
    for name, (position, radii, field) in fields.items():
        out = np.zeros((grid, grid), dtype=np.float64)
        MASL.projected_voronoi(
            out,
            position,
            field,
            radii,
            0.0,
            0.0,
            box_size,
            tracers,
            r_divisions,
            True,  # periodic
            verbose,
        )
        results[name] = out

    m_total = np.zeros((grid, grid), dtype=np.float64)
    for pos, mass, radii in [
        (masked_data.gas_position, masked_data.gas_mass, masked_data.gas_radius),
        (masked_data.dm_position, masked_data.dm_mass, masked_data.dm_radius),
        (masked_data.star_position, masked_data.star_mass, masked_data.star_radius),
        (masked_data.bh_position, masked_data.bh_mass, masked_data.bh_radius),
    ]:
        MASL.projected_voronoi(
            m_total,
            pos,
            mass,
            radii,
            0.0,
            0.0,
            box_size,
            tracers,
            r_divisions,
            True,
            verbose,
        )
    results["Mtot"] = m_total

    if masked_data.has_parts[0]:
        results["MgFe"] = results.pop("Mg") / results.pop("Fe")

        for k in ["T", "Z", "P", "Vgas"]:
            zero_mass = results["Mgas"] == 0.0
            np.divide(
                results[k],
                results["Mgas"],
                where=~zero_mass,
                out=results[k],
            )
            results[k][zero_mass] = 0.0

    if masked_data.has_parts[1] and masked_data.has_parts.sum() > 1:
        for k in ["Vcdm"]:
            zero_mass = results["Mcdm"] == 0.0
            np.divide(
                results[k],
                results["Mcdm"],
                where=~zero_mass,
                out=results[k],
            )
            results[k][zero_mass] = 0.0

    area = (box_size / grid) ** 2
    for k in ["Mgas", "Mcdm", "Mstar", "Mtot", "HI", "ne"]:
        if k in results:
            results[k] /= area

    results["z"] = data.redshift
    return results


def load_snap(snap: Path, parallelism: int) -> SnapshotData:
    """
    Assumes gadget units
    """
    with h5py.File(snap) as f:
        redshift = f["Header"].attrs["Redshift"]

        npart = f["Header"].attrs["NumPart_ThisFile"][:]
        has_parts = npart > 0

        if has_parts[0]:
            gas_position = f["PartType0/Coordinates"][:] / 1000  # Mpc / h
            gas_velocity = np.linalg.norm(
                f["PartType0/Velocities"][:], axis=1
            ) / np.sqrt(1.0 + redshift)  # km / s
            gas_mass = f["PartType0/Masses"][:] * 1e10  # Msun / h
            gas_radius = get_radii(gas_position, parallelism)
            gas_metallicity = f["PartType0/Metallicity"][:, 0] + 8e-10  # dimensionless
            gas_hI = f["PartType0/NeutralHydrogenAbundance"][:] * gas_mass  # Msun / h
            m_proton = 1.6726e-27  # kg
            Msun = 1.99e30  # kg
            kpc = 3.0857e21  # cm

            rho = f["/PartType0/Density"][:]  # (1e10 Msun/h)/(kpc/h)^3
            gas_electron = f["/PartType0/ElectronAbundance"][:]
            SFR = f["PartType0/StarFormationRate"][:]

            # formula is 0.76*ne*rho/m_proton
            # rho units are (Msun/h)/(kpc/h)^3;  1e19*2e30/(3.1e24)^3/3e-55
            factor = 1e10 * Msun / kpc**3 / m_proton

            indexes = np.where(SFR > 0.0)
            gas_electron = factor * 0.76 * gas_electron * rho  # electrons*h^2/cm^3
            gas_electron[indexes] = (
                0.0  # put electron density to 0 for star-forming particles
            )
            gas_volume = gas_mass / rho / 1e19
            gas_electron *= gas_volume

            gas_rho = f["PartType0/Density"][:]
            gas_u = f["PartType0/InternalEnergy"][:]
            gamma = 5.0 / 3.0
            gas_pressure = (gamma - 1) * gas_u * gas_rho * 1e10

            if "Temperatures" in f["PartType0"]:
                gas_temperature = f["PartType0/Temperatures"][:]
            else:
                ne = f["PartType0/ElectronAbundance"][:]
                yhelium = 0.0789
                gas_temperature = (
                    gas_u
                    * (1.0 + 4.0 * yhelium)
                    / (1.0 + yhelium + ne)
                    * 1e10
                    * (2.0 / 3.0)
                )
                BOLTZMANN = 1.38065e-16  # erg/K - NIST 2010
                PROTONMASS = 1.67262178e-24  # gram  - NIST 2010
                gas_temperature *= PROTONMASS / BOLTZMANN

            gas_magnesium = (f["PartType0/Metallicity"][:, 6] + 1e-10) * gas_mass
            gas_iron = (f["PartType0/Metallicity"][:, 10] + 1e-10) * gas_mass
        else:
            gas_position = np.empty((0, 3), dtype=np.float32)
            gas_velocity = np.empty((0, 3), dtype=np.float32)
            gas_mass = np.empty(0, dtype=np.float32)
            gas_radius = np.empty(0, dtype=np.float32)
            gas_metallicity = np.empty(0, dtype=np.float32)
            gas_hI = np.empty(0, dtype=np.float32)
            gas_electron = np.empty(0, dtype=np.float32)
            gas_temperature = np.empty(0, dtype=np.float32)
            gas_pressure = np.empty(0, dtype=np.float32)
            gas_magnesium = np.empty(0, dtype=np.float32)
            gas_iron = np.empty(0, dtype=np.float32)

        if has_parts[1]:
            dm_position = f["PartType1/Coordinates"][:] / 1000  # Mpc / h
            dm_velocity = np.linalg.norm(
                f["PartType1/Velocities"][:], axis=1
            ) / np.sqrt(1.0 + redshift)  # km / s
            dm_radius = get_radii(dm_position, parallelism)
            if "Masses" in f["PartType1"]:
                dm_mass = f["PartType1/Masses"][:] * 1e10  # Msun / h
            else:
                dm_mass = np.full_like(
                    dm_radius, f["Header"].attrs["MassTable"][1] * 1e10
                )  # Msun / h
        else:
            dm_position = np.empty((0, 3), dtype=np.float32)
            dm_velocity = np.empty((0, 3), dtype=np.float32)
            dm_mass = np.empty(0, dtype=np.float32)
            dm_radius = np.empty(0, dtype=np.float32)

        if has_parts[4]:
            star_position = f["PartType4/Coordinates"][:] / 1000  # Mpc / h
            star_mass = f["PartType4/Masses"][:] * 1e10  # Msun / h
            star_radius = np.zeros(star_position.shape[0], dtype=np.float32)
        else:
            star_position = np.empty((0, 3), dtype=np.float32)
            star_mass = np.empty(0, dtype=np.float32)
            star_radius = np.empty(0, dtype=np.float32)

        if has_parts[5]:
            bh_position = f["PartType5/Coordinates"][:] / 1000  # Mpc / h
            bh_mass = f["PartType5/Masses"][:] * 1e10  # Msun / h
            bh_radius = np.zeros(bh_position.shape[0], dtype=np.float32)
        else:
            bh_position = np.empty((0, 3), dtype=np.float32)
            bh_mass = np.empty(0, dtype=np.float32)
            bh_radius = np.empty(0, dtype=np.float32)

        box_size = f["Header"].attrs["BoxSize"] / 1000

        return SnapshotData(
            gas_position=gas_position.astype(np.float32, copy=False),
            gas_velocity=gas_velocity.astype(np.float32, copy=False),
            gas_mass=gas_mass.astype(np.float32, copy=False),
            gas_radius=gas_radius.astype(np.float32, copy=False),
            gas_metallicity=gas_metallicity.astype(np.float32, copy=False),
            gas_neutral_H=gas_hI.astype(np.float32, copy=False),
            gas_electron_density=gas_electron.astype(np.float32, copy=False),
            gas_temperature=gas_temperature.astype(np.float32, copy=False),
            gas_pressure=gas_pressure.astype(np.float32, copy=False),
            gas_magnesium=gas_magnesium.astype(np.float32, copy=False),
            gas_iron=gas_iron.astype(np.float32, copy=False),
            dm_position=dm_position.astype(np.float32, copy=False),
            dm_velocity=dm_velocity.astype(np.float32, copy=False),
            dm_mass=dm_mass.astype(np.float32, copy=False),
            dm_radius=dm_radius.astype(np.float32, copy=False),
            star_position=star_position.astype(np.float32, copy=False),
            star_mass=star_mass.astype(np.float32, copy=False),
            star_radius=star_radius.astype(np.float32, copy=False),
            bh_position=bh_position.astype(np.float32, copy=False),
            bh_mass=bh_mass.astype(np.float32, copy=False),
            bh_radius=bh_radius.astype(np.float32, copy=False),
            box_size=box_size,
            redshift=redshift,
            has_parts=has_parts,
        )


def get_radii(positions, parallelism):
    positions %= 25
    tree = KDTree(positions, boxsize=25)
    return tree.query(positions, 32 + 1, workers=parallelism)[0][:, -1]


def make_grids(data: SnapshotData, verbose: bool = False, grid: int = 256):
    result = {}
    voxel_volume = (data.box_size / grid) ** 3

    if data.has_parts[0]:
        gas_fields = {
            "T": data.gas_temperature * data.gas_mass,
            "Z": data.gas_metallicity * data.gas_mass,
            "Vgas": data.gas_velocity * data.gas_mass,
            "P": data.gas_pressure * data.gas_mass,
            "Mgas": data.gas_mass,
            "HI": data.gas_neutral_H,
            "ne": data.gas_electron_density,
            "Mg": data.gas_magnesium,
            "Fe": data.gas_iron,
        }

        field_names = []
        field_densities = []
        for name, field in gas_fields.items():
            field_density = field / (4.0 * np.pi * data.gas_radius**3 / 3.0)
            field_names.append(name)
            field_densities.append(field_density)
        densities_npy = np.stack(field_densities, axis=-1)
        with Voxelize(use_gpu=False) as v:
            output_grids = v(
                data.box_size, data.gas_position, data.gas_radius, densities_npy, grid
            )
        output_grids *= voxel_volume
        for i, name in enumerate(field_names):
            result[name] = output_grids[..., i]

        result["MgFe"] = result.pop("Mg") / result.pop("Fe")

    if data.has_parts[1]:
        dm_fields = {
            "Mcdm": data.dm_mass,
        }
        # If only dark matter
        if data.has_parts.sum() == data.has_parts[1]:
            dm_fields["Vcdm"] = data.dm_velocity * data.dm_mass

        field_names = []
        field_densities = []
        for name, field in dm_fields.items():
            field_density = field / (4.0 * np.pi * data.dm_radius**3 / 3.0)
            field_names.append(name)
            field_densities.append(field_density)
        densities_npy = np.stack(field_densities, axis=-1)
        with Voxelize(use_gpu=False) as v:
            output_grids = v(
                data.box_size, data.dm_position, data.dm_radius, densities_npy, grid
            )
        output_grids *= voxel_volume
        for i, name in enumerate(field_names):
            result[name] = output_grids[..., i]

    if data.has_parts[4]:
        m_star = np.zeros((grid, grid, grid), dtype=np.float32)
        star_position = np.copy(data.star_position)
        star_mass = np.copy(data.star_mass)
        MASL.MA(
            # data.star_position,
            star_position,
            m_star,
            data.box_size,
            "NGP",
            # W=data.star_mass,
            W=star_mass,
            verbose=verbose,
        )
        result["Mstar"] = m_star

    if data.has_parts[5]:
        m_bh = np.zeros((grid, grid, grid), dtype=np.float32)
        MASL.MA(
            data.bh_position,
            m_bh,
            data.box_size,
            "NGP",
            W=data.bh_mass,
            verbose=verbose,
        )

    # If only dark matter
    if data.has_parts.sum() == data.has_parts[1]:
        return {"Mtot": result["Mcdm"], "z": data.redshift}

    result["Mtot"] = result["Mgas"] + result["Mcdm"] + m_star + m_bh
    result["z"] = data.redshift

    return result


if __name__ == "__main__":
    main()
