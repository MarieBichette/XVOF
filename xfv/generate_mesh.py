#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Generate a 1D mesh inside a mesh.txt file
"""
import argparse
import os
from pathlib import Path
from typing import Iterable, Tuple

import numpy as np


def generate_mesh(path: Path, user_arg: Iterable[Tuple[float, int]]):
    """
    Create a mesh file

    :param path: path in with the mesh file has to be created
    :param user_arg: an array holding the length et number of points of every parts of the mesh
    """
    meshfile = path / "mesh.txt"
    with meshfile.open('w') as file_out:
        file_out.write('Mesh : initial coordinates of the nodes')
        file_out.write(os.linesep)
        file_out.write('Node Number    X coordinate [m]     Y coordinate [m]     Z coordinate [m]')
        file_out.write(os.linesep)
        file_out.write('{}         {:+10.9e}'.format(0, 0.))  # first node always at x=0
        file_out.write(os.linesep)

    total_length = 0.
    total_nb_mailles = 0
    for n_block, data in enumerate(user_arg):
        length = data[0]
        nb_mailles = data[1]
        print(f"Block number {n_block}: length = {length} - cells number = {nb_mailles}")
        bloc = np.linspace(total_length, total_length + length, int(nb_mailles) + 1)

        with meshfile.open('a') as file_out:
            for i in range(int(nb_mailles)):
                i += 1  # small hack to not write the first node which is already written
                index_node = total_nb_mailles + i
                file_out.write('{}         {:+10.9e}'.format(index_node, bloc[i]))
                file_out.write(os.linesep)

        total_length += length
        total_nb_mailles += int(nb_mailles)

    print("Generating mesh with {:} elements in repository {:}".format(total_nb_mailles, path))


def pair(arg):
    """
    Defines a pair of float, int arguments (i.e length, number of cells)
    """
    length, nb_cells = arg.split('@')
    return float(length), int(nb_cells)


def main(directory: str, meshes: Iterable[Tuple[float, int]]):
    """
    Launches the program
    """
    path = Path(directory)

    generate_mesh(path, meshes)
    print("Done !")


if __name__ == '__main__':
    # pylint: disable=invalid-name
    parser = argparse.ArgumentParser(
        description="%(prog)s generates 1D mesh files to be used by XtendedFiniteVolume",
        epilog="Example: %(prog)s base_rep 1,100 0.5,3 will generate a mesh file containing"
                "a block with a length of 1 meter and 100 cells followed by a block with a "
                "length of 50 cm and 3 cells.")
    parser.add_argument("directory", type=str,
                        help="Path toward the directory where the mesh file will be written")
    parser.add_argument("meshes", type=pair, nargs='+', metavar='length@nb_cells',
                        help="Pair length@nb_cells describint a block of mesh")
    args = parser.parse_args()

    generate_mesh(Path(args.directory), args.meshes)
    print("Done !")
