# -*- coding: utf-8 -*-
""" 
Plot error evolution in time of a field for a given item id (for a given reference field)
"""

import argparse
import pathlib
import numpy as np
import matplotlib.pyplot as plt

from xfv.post_processing.tools.hdf5_postprocessing_tools import get_field_evolution_in_time_for_item


def run():
    # ----------------------------------------------------------
    # Read instructions
    # ----------------------------------------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true", help="Increase program verbosity")
    parser.add_argument("field", help="the field to be plotted")
    parser.add_argument("item_id", type=int, help="the id of the item to look at")
    parser.add_argument("-case", action='append', nargs='+',
                        help="the path to the output repository")
    parser.add_argument("-ref", help="the path to the reference output repository")
    parser.add_argument("--output_filename", default="all_fields.hdf5",
                        help="the name of the output hdf5 band (default = all_fields.hdf5)")
    args = parser.parse_args()

    if args.case is None:
        raise ValueError("At least one case is needed. Use -case to specify cases to plot")

    if args.ref is None:
        raise ValueError("Please provide a reference case using -ref <case>")

    if args.verbose:
        print("Field : " + args.field)
        print("Item id : " + str(args.item_id))
        print("Ref : ")
        print(args.ref)
        print("Cases : ")
        print(args.case)
        print("~~~~~~~~~~~~~")

    item_id = args.item_id
    field = args.field

    # ----------------------------------------------------------
    # Prepare figure
    # ----------------------------------------------------------
    field_unit = dict()
    field_unit["Density"] = "[$kg/m^3$]"
    field_unit["Pressure"] = "[$Pa$]"
    field_unit["ArtificialViscosity"] = "[$Pa$]"
    field_unit["InternalEnergy"] = "[$J/kg$]"
    field_unit["NodeVelocity"] = "[$ms/s$]"
    field_unit["EquivalentPlasticStrainRate"] = "[$s^{-1}$]"

    plt.figure(1)
    plt.xlabel("Time [$\mu$s]")
    plt.ylabel("Diff " + field + " - ref / max(ref)")
    plt.title("Error time evolution of " + str(field) + " in item " + str(item_id), fontweight='bold')

    # ----------------------------------------------------------
    # Read reference data
    # ----------------------------------------------------------
    ref_case = args.ref
    path_to_db = pathlib.Path.cwd().joinpath("..", "tests", ref_case, args.output_filename)
    if args.verbose:
        print("Path to database : {:}".format(path_to_db))
        print("Read field " + field + " in database... ")
    # Read database :
    ref_item_history = get_field_evolution_in_time_for_item(path_to_db, item_id, field)
    ref_field = ref_item_history[:, 1]

    # ----------------------------------------------------------
    # Plot field evolution for each case
    # ----------------------------------------------------------
    for case in args.case[0]:
        if args.verbose:
            print("Case is : " + case)
        path_to_db = pathlib.Path.cwd().joinpath("..", "tests", case, args.output_filename)
        if args.verbose:
            print("Path to database : {:}".format(path_to_db))
            print("Read field " + field + " in database... ")
        # Read database :
        item_history = get_field_evolution_in_time_for_item(path_to_db, item_id, field)
        diff = np.abs(item_history[:, 1] - ref_field) / np.max(np.abs(ref_field))
        if args.verbose:
            print("Done !")
            print("~~~~~~~~~~~~~")
        # Plot field :
        plt.plot(item_history[:, 0] * 1.e+6, diff, '.-', label=case)


if __name__ == "__main__":
    run()
    # ----------------------------------------------------------
    # Show figure
    # ----------------------------------------------------------
    plt.legend(loc="best")
    plt.show()
