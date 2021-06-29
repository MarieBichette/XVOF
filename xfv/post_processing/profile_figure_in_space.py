# -*- coding: utf-8 -*-
"""
Plot the profil evolution of nodes and cell
"""

import argparse
import pathlib
import matplotlib.pyplot as plt

from xfv.post_processing.tools.hdf5_postprocessing_tools import get_field_profile_at_time
from xfv.src.output_manager.outputdatabaseexploit import OutputDatabaseExploit


def prepare_figure():
    """
    Prepare figure
    """
    plt.figure(1)
    plt.title("{:} profile".format(ARGS.field), fontweight='bold')
    plt.xlabel("Position [mm]")
    plt.ylabel("Field")


def run():
    """
    Plot field profile for each case
    """
    for case in ARGS.case[0]:
        path_to_db = pathlib.Path.cwd().joinpath(case, ARGS.output_filename)
        if ARGS.verbose:
            print("Path to database : {:}".format(path_to_db))
            print("Read field " + ARGS.field + " in database... ")
        hd_band = OutputDatabaseExploit(path_to_db)

        for current_time in ARGS.time[0]:
            # Read hdf5 band :
            coord, field_value = get_field_profile_at_time(hd_band, ARGS.field, current_time)
            if ARGS.verbose:
                print("Done !")
                print("~~~~~~~~~~~~~")
            # Plot field :
            plt.plot(coord * 1.e+03, field_value, label=case)

        if (ARGS.write_data):
            data_path="{:s}Profile_{:s}_{:3.1e}.dat".format(case, ARGS.field, current_time)
            with open(data_path, "w") as file_object:
                for x_data, y_data in zip(coord, field_value):
                    file_object.write("{:20.18g}\t{:20.18g}\n".format(x_data, y_data))
            print("Data written in {:s}".format(data_path))


if __name__ == "__main__":
    # -----------------------------------------
    # Read instructions
    # -----------------------------------------
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("-v", "--verbose", action="store_true", help="Increase program verbosity")
    PARSER.add_argument("field", help="the field to be plotted")
    PARSER.add_argument("time", type=float, action='append', nargs='+',
                        help="the list of times to look at")
    PARSER.add_argument("-case", action='append', nargs='+',
                        help="the path to the output repository")
    PARSER.add_argument("--output_filename", default="all_fields.hdf5",
                        help="the name of the output hdf5 band (default = all_fields.hdf5)")
    PARSER.add_argument("--write_data", action="store_true", help="To write data in an output file")
    ARGS = PARSER.parse_args()

    if ARGS.case is None:
        raise ValueError("At least one case is needed. Use -case to specify cases to plot")

    if ARGS.time is None:
        raise ValueError("At least one time is needed. Use -time to specify times to plot")

    if ARGS.verbose:
        print("Cases : ")
        print(ARGS.case)
        print("~~~~~~~~~~~~~")

    # -----------------------------------------
    # Figure preparation
    # -----------------------------------------
    prepare_figure()

    # -----------------------------------------
    # Run post processing
    # -----------------------------------------
    run()

    # -----------------------------------------
    # Show figure
    # -----------------------------------------
    plt.legend(loc="best")
    plt.show()
