#!/bin/env python
import math
import sys

import h5py
import numpy as np
import pygadgetreader


def main():
    if len(sys.argv) != 3:
        print(f"USAGE: {sys.argv[0]} [source file] [destination file]")
        return

    source = sys.argv[1]
    destination = sys.argv[2]

    copy_ic(source, destination)


def copy_ic(source: str, destination: str):
    with h5py.File(destination, "w") as f:
        header = {}
        original_header = pygadgetreader.readheader(source, "header")
        h = original_header["h"]
        a = original_header["time"]
        header["BoxSize"] = original_header["boxsize"]
        header["Dimension"] = 3
        header["Flag_Entropy_ICs"] = 0
        header["MassTable"] = [0.0] * 7
        header["NumFilesPerSnapshot"] = 1
        header["NumPart_ThisFile"] = [
            original_header["ngas"],
            original_header["ndm"],
        ] + [0] * 5
        header["NumPart_Total"] = [
            original_header["ngas"] & 0xFFFFFFFF,
            original_header["ndm"] & 0xFFFFFFFF,
        ] + [0] * 5
        header["NumPart_Total_HighWord"] = [
            original_header["ngas"] >> 32,
            original_header["ndm"] >> 32,
        ] + [0] * 5
        header["Redshift"] = original_header["redshift"]
        header["Time"] = original_header["time"]

        f_header = f.create_group("Header")
        f_header.attrs.update(header)

        fields = []
        if original_header["ndm"] > 0:
            fields += [
                ("pos", "dm", "PartType1/Coordinates"),
                ("vel", "dm", "PartType1/Velocities"),
                ("mass", "dm", "PartType1/Masses"),
                ("pid", "dm", "PartType1/ParticleIDs"),
            ]
        if original_header["ngas"] > 0:
            fields += [
                ("pos", "gas", "PartType0/Coordinates"),
                ("vel", "gas", "PartType0/Velocities"),
                ("mass", "gas", "PartType0/Masses"),
                ("pid", "gas", "PartType0/ParticleIDs"),
            ]

        for gadget_name, gadget_type, swift_name in fields:
            data = pygadgetreader.readsnap(source, gadget_name, gadget_type)
            if gadget_name == "pos":
                # Make sure coordinates are withing the size of the box
                # Wrap them around if they're not.
                data %= header["BoxSize"]
            f[swift_name] = data

        # Set smoothing lengths to mean inter-particle separation
        hsm = header["BoxSize"] / np.cbrt(original_header["ngas"])
        f["PartType0/SmoothingLength"] = hsm * np.ones(original_header["ngas"])

        # Write internal energies
        # based on code from monofonic
        gamma = 5 / 3
        YHe = 0.2457519853817943
        M_gas = f["PartType0/Masses"][:].sum()
        M_dm = f["PartType1/Masses"][:].sum()
        Omega_b = M_gas / (M_gas + M_dm) * original_header["O0"]
        Tcmb0 = 2.7255

        npol = 1 / (gamma - 1) if math.fabs(1.0 - gamma) > 1e-7 else 1
        unitv = 1e5
        adec = 1.0 / (160.0 * (Omega_b * h * h / 0.022) ** (2.0 / 5.0))
        Tini = Tcmb0 / a if a < adec else Tcmb0 / a / a * adec
        mu = 4.0 / (8.0 - 5.0 * YHe) if Tini > 1e4 else 4.0 / (1.0 + 3.0 * (1.0 - YHe))
        energy = 1.3806e-16 / 1.6726e-24 * Tini * npol / mu / unitv / unitv
        energies = np.full(original_header["ngas"], energy)
        f["PartType0/InternalEnergy"] = energies

        # Make sure all particles are in the box.
        # Sometimes the conversion results in particles just barely outside the box.
        f["PartType0/Coordinates"][:] %= header["BoxSize"]
        f["PartType1/Coordinates"][:] %= header["BoxSize"]


if __name__ == "__main__":
    main()
