import h5py
import argparse
from pathlib import Path
from astropy.table import Table

import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("ndskel", type=Path)
    parser.add_argument("subfind", type=Path)
    parser.add_argument("param_file", type=Path, nargs="?", default=None)
    parser.add_argument("--target", type=Path, default=Path.cwd())
    parser.add_argument("--boxsize", type=float, default=25000)
    parser.add_argument("--grid", type=int, default=256)

    args = parser.parse_args()
    ndskel: Path = args.ndskel.resolve()
    subfind: Path = args.subfind.resolve()
    param_file: Path | None = args.param_file
    if param_file is not None:
        param_file = param_file.resolve()
    target: Path = args.target
    boxsize: float = args.boxsize
    grid: int = args.grid

    with h5py.File(subfind) as f:
        subhalo_positions: np.ndarray = f["Subhalo/SubhaloPos"][:]

    multiply_factor = boxsize / grid
    skeleton_data = parse_ndskl_ascii(ndskel)
    critical_points = (
        np.array([cp["position"] for cp in skeleton_data["critical_points"]["data"]])
        * multiply_factor
    )

    closest_halo_indices, closest_distances = find_closest_halo(
        critical_points, subhalo_positions
    )

    filaments_df = extract_filaments_data(skeleton_data, multiply_factor, boxsize)
    critical_points_df = extract_critical_points_data(
        skeleton_data, closest_halo_indices, closest_distances, multiply_factor
    )

    structured_data = {
        "critical_points": table2dict(critical_points_df),
        "filaments": table2dict(filaments_df),
    }

    if param_file is not None:
        sim_dir: Path = ndskel.relative_to(param_file.parent).parents[-2]
        cosmoastrotable: Table = Table.read(param_file, format="ascii.basic")
        name = cosmoastrotable.colnames[0]
        cosmoastrotable.add_index(name)
        row = cosmoastrotable.loc[str(sim_dir)]
        row_dict = {k: row[k] for k in row.colnames if k != name and k != "seed"}
        structured_data["header"] = {}
        structured_data["header"]["parameters"] = row_dict

    save_to_hdf5_disperse(
        structured_data,
        target / (ndskel.stem.split(".", maxsplit=1)[0] + "-output_file_090.hdf5"),
    )


def table2dict(table: Table):
    result = {}
    for col in table.colnames:
        column = table[col]
        dtype = column.dtype
        if dtype in (np.float64, np.int64):
            result[col] = np.array(column)
        elif dtype == np.dtype("O"):
            result[col] = list(column)
        else:
            print(col, dtype)
    return result


def find_closest_halo(critical_points, halo_positions):
    """
    Finds the closest halo to each critical point in 3D space.

    Parameters
    ----------
    critical_points : np.ndarray
        Array of critical point positions, shape (N, 3) for X, Y, Z.
    halo_positions : np.ndarray
        Array of halo positions, shape (M, 3) for X, Y, Z.

    Returns
    -------
    closest_halo_indices : np.ndarray
        Array of indices of the closest halo to each critical point.
    closest_distances : np.ndarray
        Array of distances to the closest halo for each critical point.
    """
    closest_halo_indices = np.zeros(len(critical_points), dtype=int)
    closest_distances = np.zeros(len(critical_points), dtype=float)

    for i, cp in enumerate(critical_points):
        # Compute the distance to all haloes
        delta = halo_positions - cp
        distances = np.linalg.norm(delta, axis=1)

        # Find the index of the minimum distance
        closest_idx = np.argmin(distances)
        closest_halo_indices[i] = closest_idx
        closest_distances[i] = distances[closest_idx]

    return closest_halo_indices, closest_distances


def parse_ndskl_ascii(fname: Path):
    """
    Parses an NDskl ASCII file into structured data.

    Parameters
    ----------
    filename : str
        The path to the NDskl ASCII file.

    Returns
    -------
    data : dict
        A dictionary containing structured data from the NDskl ASCII file.
    """
    with open(fname, "r") as f:
        lines = f.readlines()

    data = {
        "header": {},
        "critical_points": {"ncrit": 0, "data": []},
        "filaments": {"nfil": 0, "data": []},
        "critical_points_data": {"NF": 0, "fields": [], "data": []},
        "filaments_data": {"NF": 0, "fields": [], "data": []},
    }
    section = None
    ndims = None
    i = 0

    while i < len(lines):
        line = lines[i].strip()
        if not line or line.startswith("#"):
            i += 1
            continue

        if line.startswith("ANDSKEL"):
            section = "header"
            i += 1
            continue

        if line == "[CRITICAL POINTS]":
            section = "critical_points"
            i += 1
            data["critical_points"]["ncrit"] = int(lines[i].strip())
            i += 1
            continue

        if line == "[FILAMENTS]":
            section = "filaments"
            i += 1
            data["filaments"]["nfil"] = int(lines[i].strip())
            i += 1
            continue

        if line == "[CRITICAL POINTS DATA]":
            section = "critical_points_data"
            i += 1
            data["critical_points_data"]["NF"] = int(lines[i].strip())
            for _ in range(data["critical_points_data"]["NF"]):
                i += 1
                data["critical_points_data"]["fields"].append(lines[i].strip())
            i += 1
            continue

        if line == "[FILAMENTS DATA]":
            section = "filaments_data"
            i += 1
            data["filaments_data"]["NF"] = int(lines[i].strip())
            for _ in range(data["filaments_data"]["NF"]):
                i += 1
                data["filaments_data"]["fields"].append(lines[i].strip())
            i += 1
            continue

        if section == "header":
            key, value = parse_header_line(line)
            if key:
                data["header"][key] = value
                if key == "ndims":
                    ndims = value

        elif section == "critical_points":
            cp_data = line.split()
            cp = {
                "type": int(cp_data[0]),
                "position": np.array(
                    list(map(float, cp_data[1 : 1 + ndims]))[::-1]
                ),  # Reverse to Z, Y, X
                "value": float(cp_data[1 + ndims]),
                "pairID": int(cp_data[2 + ndims]),
                "boundary": int(cp_data[3 + ndims]),
                "filaments": [],
            }
            i += 1
            nfil = int(lines[i].strip())
            for _ in range(nfil):
                i += 1
                destId, filId = map(int, lines[i].strip().split())
                cp["filaments"].append({"destId": destId, "filId": filId})
            data["critical_points"]["data"].append(cp)

        elif section == "filaments":
            fil_data = line.split()
            filament = {
                "CP1": int(fil_data[0]),
                "CP2": int(fil_data[1]),
                "nSamp": int(fil_data[2]),
                "sampling_points": [],
            }
            for _ in range(filament["nSamp"]):
                i += 1
                samp_point = np.array(
                    list(map(float, lines[i].strip().split()))[::-1]
                )  # Reverse to Z, Y, X
                filament["sampling_points"].append(samp_point)
            data["filaments"]["data"].append(filament)

        elif section == "critical_points_data":
            cp_data_values = list(map(float, line.split()))
            data["critical_points_data"]["data"].append(cp_data_values)

        elif section == "filaments_data":
            fil_data_values = list(map(float, line.split()))
            data["filaments_data"]["data"].append(fil_data_values)

        i += 1

    return data


def parse_header_line(line):
    """
    Parses a line from the header section.

    Parameters
    ----------
    line : str
        The line from the header to parse.

    Returns
    -------
    key : str
        The key indicating the type of data parsed (e.g., 'ndims', 'bbox').
    value : int, list
        The parsed value associated with the key.
    """
    parts = line.split()
    if len(parts) == 1 and parts[0].isdigit():
        return "ndims", int(parts[0])
    elif line.startswith("BBOX"):
        bbox_str = (
            line.replace("BBOX", "").replace("[", "").replace("]", "").replace(",", " ")
        )
        bbox_parts = list(map(float, bbox_str.split()))
        ndims = (
            len(bbox_parts) // 2
        )  # Assume ndims is half of the parts, since there are two sets of ndims values
        bbox_origin = bbox_parts[:ndims]
        bbox_delta = bbox_parts[ndims:]
        return "bbox", [bbox_origin, bbox_delta]
    return None, None


def extract_filaments_data(data, multiply_factor, boxsize):
    """
    Extracts the properties and positions of all filaments from the parsed NDskl data.

    Parameters
    ----------
    data : dict
        The parsed NDskl data containing filaments and other structures.
    multiply_factor : float
        Factor to multiply positions by to convert to the correct physical units.

    Returns
    -------
    filaments_df : pd.DataFrame
        A pandas DataFrame containing the properties and positions of all filaments.
    """
    filaments_list = []

    for filament, filament_data in zip(
        data["filaments"]["data"], data["filaments_data"]["data"]
    ):
        # Apply multiply_factor to sampling points
        sampling_points = np.array(filament["sampling_points"]) * multiply_factor

        filament_info = {
            "CP1_index": filament.get("CP1", None),  # Index of the first critical point
            "CP2_index": filament.get(
                "CP2", None
            ),  # Index of the second critical point
            "nSamp": filament.get("nSamp", None),  # Number of sampling points
            "sampling_points": sampling_points,  # Scaled list of sampling points (positions)
            "length": compute_filament_length(
                sampling_points, boxsize
            ),  # Length of the filament
            "field_value": filament_data[0],
            "orientation": filament_data[1],
            "cell": filament_data[2],
            "log_field_value": filament_data[3],
            "type": filament_data[4],
        }
        filaments_list.append(filament_info)

    return Table(filaments_list)


def compute_filament_length(sampling_points, box_size):
    """
    Computes the length of a filament accounting for periodic boundary conditions.

    Parameters
    ----------
    sampling_points : np.ndarray
        Array of filament points, shape (N, 3) for X, Y, Z coordinates.
    box_size : float
        The size of the periodic box.

    Returns
    -------
    length : float
        The total length of the filament.
    """
    length = 0.0

    # Iterate over pairs of consecutive points
    for i in range(len(sampling_points) - 1):
        p1 = sampling_points[i]
        p2 = sampling_points[i + 1]

        # Compute the distance vector between p1 and p2
        delta = p2 - p1

        # Apply periodic boundary conditions
        delta = delta - np.round(delta / box_size) * box_size

        # Compute the Euclidean distance
        segment_length = np.linalg.norm(delta)

        # Sum the segment length to the total length
        length += segment_length

    return length


def extract_critical_points_data(
    data, closest_halo_indices, closest_distances, multiply_factor
):
    """
    Extracts the positions and properties of all critical points from the parsed NDskl data, including the closest halo information.

    Parameters
    ----------
    data : dict
        The parsed NDskl data containing critical points and other structures.
    closest_halo_indices : np.ndarray
        Array of indices of the closest halo to each critical point.
    closest_distances : np.ndarray
        Array of distances to the closest halo for each critical point.
    multiply_factor : float
        Factor to multiply positions by to convert to the correct physical units.

    Returns
    -------
    critical_points_df : pd.DataFrame
        A pandas DataFrame containing the positions and properties of all critical points.
    """
    critical_points_list = []

    for i, (cp, cpd) in enumerate(
        zip(data["critical_points"]["data"], data["critical_points_data"]["data"])
    ):
        filaments = cp.get(
            "filaments", []
        )  # Get the list of filaments or an empty list if missing
        nfil = len(filaments)  # Number of filaments is the length of this list

        # Apply multiply_factor to position coordinates
        position = np.array(cp["position"]) * multiply_factor

        cp_info = {
            "type": cp.get("type", None),
            "position_x": position[0],
            "position_y": position[1],
            "position_z": position[2],
            "value": cp.get("value", None),
            "pairID": cp.get("pairID", None),
            "boundary": cp.get("boundary", None),
            "nfil": nfil,  # Number of filaments
            "filaments": filaments,  # List of filaments
            "closest_halo_index": closest_halo_indices[i],
            "closest_halo_distance": closest_distances[i],
            "persistence_ratio": cpd[0],
            "persistence_nsigmas": cpd[1],
            "persistence": cpd[2],
            "persistence_pair": cpd[3],
            "parent_index": cpd[4],
            "parent_log_index": cpd[5],
            "log_field_value": cpd[6],
            "field_value": cpd[7],
            "cell": cpd[8],
        }
        critical_points_list.append(cp_info)

    # Convert the list of dictionaries to a pandas DataFrame
    return Table(critical_points_list)


def save_to_hdf5_disperse(data, hdf5_filename):
    """
    Saves the parsed NDskl data to an HDF5 file.

    Parameters
    ----------
    data : dict
        Parsed NDskl data.
    hdf5_filename : str
        The path to the output HDF5 file.
    """
    # Analyze the dictionary structure
    print("Analyzing dictionary structure...")
    analyze_dict_structure(data)

    with h5py.File(hdf5_filename, "w") as hf:

        def recursively_add_dict_to_group(group, dictionary):
            """
            Recursively saves a dictionary to an HDF5 group.

            Parameters
            ----------
            group : h5py.Group
                HDF5 group to save the data into.
            dictionary : dict
                Dictionary containing the data to save.
            """
            for key, value in dictionary.items():
                # Check for lists longer than 10
                if isinstance(value, list) and len(value) > 10:
                    print(
                        f"Processing key: {key}, type: {type(value)}, list is too long to print..."
                    )
                # Check for arrays longer than 10 rows
                elif isinstance(value, np.ndarray) and value.shape[0] > 10:
                    print(
                        f"Processing key: {key}, type: {type(value)}, shape: {value.shape}, array is too long to print..."
                    )
                # Check for dictionaries
                elif isinstance(value, dict):
                    print(
                        f"Processing key: {key}, type: {type(value)} (dictionary with keys: {list(value.keys())})"
                    )
                # Default case for other types
                else:
                    print(f"Processing key: {key}, type: {type(value)}, value: {value}")

                if isinstance(value, dict):
                    # If the value is a dictionary, create a subgroup
                    print(f"Creating subgroup for key: {key}")
                    subgroup = group.create_group(key)
                    recursively_add_dict_to_group(subgroup, value)
                elif isinstance(value, list):
                    try:
                        if all(isinstance(i, (int, float, str)) for i in value):
                            print(f"Saving list to dataset: {key}")
                            group.create_dataset(
                                key,
                                data=np.array(
                                    value,
                                    dtype=np.string_
                                    if isinstance(value[0], str)
                                    else np.float64,
                                ),
                            )
                        elif all(isinstance(i, dict) for i in value):
                            print(f"Handling list of dictionaries: {key}")
                            subgroup = group.create_group(key)
                            for idx, item in enumerate(value):
                                item_group = subgroup.create_group(f"item_{idx}")
                                recursively_add_dict_to_group(item_group, item)
                        elif all(
                            isinstance(i, list) and all(isinstance(j, dict) for j in i)
                            for i in value
                        ):
                            print(f"Handling list of lists of dictionaries: {key}")
                            subgroup = group.create_group(key)
                            for idx, sublist in enumerate(value):
                                sublist_group = subgroup.create_group(
                                    f"critical_point_{idx}"
                                )
                                for jdx, subdict in enumerate(sublist):
                                    subdict_group = sublist_group.create_group(
                                        f"filament_{jdx}"
                                    )
                                    recursively_add_dict_to_group(
                                        subdict_group, subdict
                                    )
                        else:
                            print(f"Handling complex list: {key}")
                            subgroup = group.create_group(key)
                            for idx, item in enumerate(value):
                                if isinstance(item, (int, float, str, np.ndarray)):
                                    subgroup.create_dataset(
                                        f"filament_{idx}", data=item
                                    )
                                else:
                                    subgroup.attrs[f"filament_{idx}"] = str(item)
                    except Exception as e:
                        print(f"Error saving list for key {key}: {e}")
                elif isinstance(value, (int, float, np.ndarray)):
                    # Save numerical data directly
                    print(f"Saving numeric data: {key}")
                    group.create_dataset(key, data=value)
                elif isinstance(value, str):
                    # Save strings as attributes
                    print(f"Saving string as attribute: {key}")
                    group.attrs[key] = value
                else:
                    # Serialize unsupported types to strings
                    print(f"Serializing unsupported type for key {key}: {type(value)}")
                    group.attrs[key] = str(value)

        # Start the recursive saving process
        recursively_add_dict_to_group(hf, data)

    print(f"Data successfully saved to {hdf5_filename}")


def analyze_dict_structure(data, level=0, parent_key="root"):
    """
    Recursively analyze the structure of a dictionary, logging the type of each value.

    Parameters
    ----------
    data : dict
        Dictionary to analyze.
    level : int
        Current depth in the dictionary hierarchy (used for indentation).
    parent_key : str
        Parent key for the current dictionary (for clarity in output).
    """
    indent = "  " * level
    for key, value in data.items():
        full_key = f"{parent_key}/{key}"
        if isinstance(value, dict):
            print(f"{indent}{full_key} (dict): {len(value)} keys")
            analyze_dict_structure(value, level + 1, full_key)
        elif isinstance(value, list):
            print(f"{indent}{full_key} (list): {len(value)} items")
            if len(value) > 0:
                item_types = {type(item).__name__ for item in value}
                print(f"{indent}  List contains types: {item_types}")
        else:
            print(f"{indent}{full_key} ({type(value).__name__}): {value}")


if __name__ == "__main__":
    main()
