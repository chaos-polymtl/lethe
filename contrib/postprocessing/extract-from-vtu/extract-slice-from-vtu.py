"""
SPDX-FileCopyrightText: Copyright (c) 2019-2025 The Lethe Authors
SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

This script extracts a slice from a .vtu file or all the .vtu files linked to a .pvd file. The slice is defined by a point and a normal direction and is saved in a new .vtu file.

To use: python3 extract-slice-from-vtu.py <file> --origin <x>,<y>,<z> --normal <u>,<v>,<w> --np <number of processes>
"""

"""
IMPORTS
"""

import sys
import numpy as np
import pyvista as pv
import argparse
import os
import xml.etree.ElementTree as ET
from multiprocessing import Pool
from functools import partial
from tqdm import tqdm

"""
FUNCTIONS
"""

def get_pvtus_path_from_pvd(pvd_file_path: str) -> list:
    """
    Extract the .vtu files if a .pvd file is provided.

    :param pvd_file_path: the path to the .pvd file
    :return: a list of the paths to the .vtu files
    """

    # Extract .pvtu files from the .pvd file
    tree = ET.parse(pvd_file_path)
    root = tree.getroot()
    pvtu_files = [dataset.get("file") for dataset in root.findall(".//DataSet")]

    # Extract .pvtu file
    pvtu_files_complete = []
    for pvtu_file in pvtu_files:
        pvtu_files_complete.append(os.path.join(os.path.dirname(pvd_file_path), pvtu_file))

    return pvtu_files_complete

def extract_slice(file_path: list, origin: np.array, normal: np.array, save_path: str) -> None:
    """
    Create a slice from a .vtu file using the pyvista library and save it in a new .vtk file.

    :param file_path: the path to the .vtu or .pvtu file
    :param origin: the coordinates of a point on the slice
    :param normal: the normal direction to the slice
    :param save_path: the path to the sliced vtu save folder
    """

    # slice the file
    sliced_solution = pv.read(file_path).slice(normal=normal, origin=origin, generate_triangles=True)
    sliced_solution = sliced_solution.cast_to_unstructured_grid()

    # If the slice is empty (the slice does not intersect the domain), return
    if sliced_solution.n_points == 0:
        print(f"Slice is empty for file {file_path}")
        return

    # Modify the file name and save the file
    file_name = os.path.basename(file_path)

    # If a .pvtu file was provided, change the extension to .vtu, because this tool writes the solution in .vtu file format
    file_name=file_name.replace(".pvtu", ".vtu")

    sliced_solution.save(save_path+"/sliced-" + file_name)

    return

"""
MAIN
"""

def run()-> None:
    """
    Run the script.
    """

    parser = argparse.ArgumentParser(description="Extract a slice from a vtu file")
    parser.add_argument("file", help="Path to a .vtu or a .pvd file")
    parser.add_argument("--origin", help="The coordinates of a point on the slice", default="0,0,0")
    parser.add_argument("--normal", help="The normal direction to the slide", default="0,1,0")
    parser.add_argument("--np", type=int, help="The number of processes", default="1")
    parser.add_argument("--save_folder", type=str, help="Path to where the sliced .vtu are saved", default="./")
    args = parser.parse_args()

    # set the origin and normal
    origin = np.array([float(x) for x in args.origin.split(",")])
    normal = np.array([float(x) for x in args.normal.split(",")])

    print("Slicing vtus using point: ", origin, " and normal: ", normal)

    if (args.file.endswith(".vtu") or args.file.endswith(".pvtu")):
        extract_slice(args.file, origin=origin, normal=normal, save_path=args.save_folder)
    elif args.file.endswith(".pvd"):
        pvtu_files_path = get_pvtus_path_from_pvd(args.file)
        if args.np > 1:
            pool = Pool(args.np)
            list(tqdm(pool.imap(partial(extract_slice, origin=origin, normal=normal, save_path=args.save_folder), pvtu_files_path), 
                        total=len(pvtu_files_path), desc="Processing PVTU files"))
        else:
            for pvtu_file_path in tqdm(pvtu_files_path, desc="Processing PVTU files"):
                extract_slice(pvtu_file_path, origin, normal, args.save_folder)
        
        # Write the sliced .pvd
        sliced_pvd_file_name = args.save_folder+"/sliced-"+ os.path.basename(args.file)
        sliced_pvd_file = open(sliced_pvd_file_name, "w")
        pvd_file = open(args.file,"r")

        pvd_file_lines = pvd_file.read()
        pvd_file_lines = pvd_file_lines.replace('file="','file="sliced-')
        pvd_file_lines = pvd_file_lines.replace(".pvtu", ".vtu")
        sliced_pvd_file.write(pvd_file_lines)

        sliced_pvd_file.close()
        pvd_file.close()

    else:
        print("Invalid file format, please provide a .vtu or a .pvd file")

if __name__ == "__main__":
    run()
