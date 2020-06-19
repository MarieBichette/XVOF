# -*- coding: utf-8 -*-
"""
Plot time evolution of a field for a given item id
"""

import argparse
import pathlib
import matplotlib.pyplot as plt

from xfv.post_processing.tools.hdf5_postprocessing_tools import get_field_evolution_in_time_for_item


def run():
    """
    Run the postprocessing program
    """
    # ----------------------------------------------------------
    # Read instructions
    # ----------------------------------------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true", help="Increase program verbosity")
    parser.add_argument("field", help="the field to be plotted")
    parser.add_argument("item_id", type=int, help="the id of the item to look at")
    parser.add_argument("-case", action='append', nargs='+',
                        help="the path to the output repository")
    parser.add_argument("--output_filename", default="all_fields.hdf5",
                        help="the name of the output hdf5 band (default = all_fields.hdf5)")
    parser.add_argument("--write_data", action="store_true", help="To write data in an output file")
    args = parser.parse_args()

    if args.case is None:
        raise ValueError("At least one case is needed. Use -case to specify cases to plot")

    if args.verbose:
        print("Field : " + args.field)
        print("Item id : " + str(args.item_id))
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
    field_unit["Porosity"] = "[$.$]"
    field_unit["ShearModulus"] = "[$Pa$]"
    field_unit["YieldStress"] = "[$Pa$]"
    
    plt.figure(1)
    plt.xlabel("Time [mus]")
    plt.ylabel(field + field_unit[field])
    plt.title("Time evolution of " + str(field) + " in item " + str(item_id), fontweight='bold')

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
        if args.verbose:
            print("Done !")
            print("~~~~~~~~~~~~~")
        # Plot field :
        plt.plot(item_history[:, 0] * 1.e+6, item_history[:, 1], '.-', label=case)
        if(args.write_data):
            data_path="{:s}Field_evolution_{:s}_{:d}.dat".format(case,field,item_id)
            with open(data_path, "w") as file_object:
                for x_data, y_data in zip(item_history[:, 0], item_history[:, 1]):
                    file_object.write("{:20.18g}\t{:20.18g}\n".format(x_data,y_data))
            print("Data written in {:s}".format(data_path))


if __name__ == "__main__":
    run()
    # ----------------------------------------------------------
    # Show figure
    # ----------------------------------------------------------
    plt.legend(loc="best")
    plt.show()
