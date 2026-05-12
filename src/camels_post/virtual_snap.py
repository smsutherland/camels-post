import argparse
import typing as T
from pathlib import Path

import h5py
import hdf5plugin as _  # noqa
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "snaps",
        help="First snapshot in multi-file snapshot or list of all snapshots.",
        type=Path,
        nargs="+",
    )
    parser.add_argument("destination", help="Destination virtual snapshot", type=Path)
    parser.add_argument(
        "--force",
        "-f",
        help="Override if destination already exists",
        action="store_true",
    )
    parser.add_argument(
        "--resolve",
        "-r",
        help="Resolve snapshot names to absolute paths",
        action="store_true",
    )
    args = parser.parse_args()
    snaps: list[Path] = args.snaps
    dst: Path = args.destination
    force: bool = args.force
    resolve: bool = args.resolve

    if resolve:
        snaps = [s.resolve() for s in snaps]
    else:
        snaps = [
            s.relative_to(dst.parent) if s.is_relative_to(dst.parent) else s.resolve()
            for s in snaps
        ]

    if len(snaps) == 1:
        if not str(snaps[0]).endswith(".0.hdf5"):
            print(f"ERROR: Unrecognized snapshot name format {snaps[0]}")
            print("Please provide all snapshot names directly as arguments")
            exit(1)
        i = 1
        base = snaps[0].parent
        base_name = snaps[0].stem[:-2]
        while (f := base / (base_name + f".{i}.hdf5")).exists():
            i += 1
            snaps.append(f)

    non_existent = [s for s in snaps if not s.exists()]
    if len(non_existent) > 0:
        print("ERROR: Snapshots")
        for s in non_existent:
            print(" ", s)
        print("do not exist!")
        exit(1)

    if dst.exists():
        if force:
            dst.unlink()
        else:
            print("ERROR:", dst, "already exists!")
            exit(1)

    print("combining")
    for s in snaps:
        print(" ", s)
    print("into", dst)

    combine(snaps, dst)


def combine(snaps: list[Path], dst: Path):
    main_snap = snaps[0]
    with h5py.File(main_snap) as f:
        header = dict(f["Header"].attrs)
    header["Virtual"] = 1
    header["NumFilesPerSnapshot"] = 1
    num_part_total = header["NumPart_Total"]
    if np.any(header["NumPart_Total_HighWord"] > 0):
        print("Warning: NumPart > 2^32. Switching to 64 bit NumPart!")
        num_part_dtype = np.uint64
        num_part_total = num_part_total.astype(np.uint64)
        num_part_total += header["NumPart_Total_HighWord"].astype(np.uint64) << 32
    else:
        num_part_dtype = np.uint32
    header["NumPart_ThisFile"] = num_part_total

    ptypes = np.flatnonzero(num_part_total)

    datasets: dict[tuple[int, str], h5py.VirtualLayout] = {}
    with h5py.File(snaps[0]) as f:
        for p in ptypes:
            g = f[f"PartType{p}"]
            for k in g.keys():
                ds = g[k]
                datasets[(p, k)] = h5py.VirtualLayout(
                    shape=(num_part_total[p],) + ds.shape[1:],
                    dtype=ds.dtype,
                )

    num_part = np.zeros_like(num_part_total, dtype=num_part_dtype)
    for s in snaps:
        with h5py.File(s) as hs:
            npart_this_file = (
                hs["Header"].attrs["NumPart_ThisFile"][:].astype(num_part_dtype)
            )
            for (p, ds_name), vds in datasets.items():
                vds[num_part[p] : num_part[p] + npart_this_file[p]] = (
                    h5py.VirtualSource(hs[f"PartType{p}"][ds_name])
                )
            num_part += npart_this_file

    with h5py.File(dst, "w") as dstf:
        dstf: h5py.File = dstf

        file_header: h5py.Group = dstf.create_group("Header")
        file_header.attrs.update(header)

        for p in ptypes:
            dstf.create_group(f"PartType{p}")

        for (p, ds_name), vds in datasets.items():
            g: h5py.Group = dstf[f"PartType{p}"]
            g.create_virtual_dataset(ds_name, vds)
