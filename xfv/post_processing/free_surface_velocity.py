# -*- coding: utf-8 -*-
"""
Plot the free surface velocity eventually with experimental data
"""

import argparse
import pathlib
import numpy as np
import matplotlib.pyplot as plt

from xfv.post_processing.tools.hdf5_postprocessing_tools import get_field_evolution_in_time_for_item


def run():
    """
    Run post processing program
    """
    # ----------------------------------------------------------
    # Read instructions
    # ----------------------------------------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true", help="Increase program verbosity")
    parser.add_argument("-case", action='append', nargs='+',
                        help="the path to the output repository from /xfv/tests/")
    parser.add_argument("-experimental_data", help="the path to experimental data to plot")
    parser.add_argument("--output_filename", default="all_fields.hdf5",
                        help="the name of the output hdf5 band (default = all_fields.hdf5)")
    parser.add_argument("--shift_t0", action="store_true",
                        help="Shift time origin to put t0 when the velocity signal arrives on "
                             "the free surface velocity")
    parser.add_argument("--write_data", default="free_surface_velocity.dat",
                        help="Name of the output file")
    args = parser.parse_args()

    if args.case is None:
        raise ValueError("At least one case is needed. Use -case to specify cases to plot")

    if args.verbose:
        print("Cases : ")
        print(args.case)
        if args.experimental_data is not None:
            print("Experimental data : " + args.experimental_data)
        print("~~~~~~~~~~~~~")

    exp_data = args.experimental_data

    # ----------------------------------------------------------
    # Prepare figure
    # ----------------------------------------------------------
    plt.figure(1)
    plt.title("Evolution of the free surface velocity", fontweight='bold', fontsize=18)
    plt.xlabel("Time [mus]", fontsize=16)
    plt.ylabel("Free surface velocity [m/s]", fontsize=16)

    # ----------------------------------------------------------
    # Plot free surface velocity for each case
    # ----------------------------------------------------------
    for case in args.case[0]:
        if args.verbose:
            print("Case is : " + case)
        path_to_db = pathlib.Path.cwd().joinpath(case, args.output_filename)
        if args.verbose:
            print("Path to database : {:}".format(path_to_db))
            print("Read VelocityField in database... ")
        # Read database :
        # Free surface is the last node => index -1 in Numpy array
        item_history = get_field_evolution_in_time_for_item(path_to_db, -1, "ClassicalNodeVelocity")
        # NodeVelocity => calls the true_field_extraction with NodeCoordinates reading
        # ClassicalNodeVelocity => calls field_extraction to get the field without NodeCoordinates extraction
        if args.verbose:
            print("Done !")
            print("~~~~~~~~~~~~~")
        # Plot velocity with t0 = detection of free surface movement:
        time = item_history[:, 0]
        velocity = item_history[:, 1]
        time_0 = 0.
        if args.shift_t0:
            time_0 = time[velocity > 1][0]  # 1st time where non zero velocity
            if args.verbose:
                print("New t0 is : " + str(time_0))
        plt.plot((time - time_0) * 1.e+6, velocity, label=case)
        if (args.write_data):
            data_path=case + args.write_data
            with open(data_path, "w") as file_object:
                for x_data, y_data in zip(time, velocity):
                    file_object.write("{:20.18g}\t{:20.18g}\n".format(x_data, y_data))
            print("Data written in {:s}".format(data_path))

    # ----------------------------------------------------------
    # Plot experimental data
    # ----------------------------------------------------------
    if exp_data is not None:
        experimental_velocity = np.loadtxt(exp_data)
        plt.plot(experimental_velocity[:, 0], experimental_velocity[:, 1], "--",
                 color="black", label="Experiment")


if __name__ == "__main__":
    run()
    # ----------------------------------------------------------
    # Show figure
    # ----------------------------------------------------------
    plt.legend(loc="best")
    plt.show()
