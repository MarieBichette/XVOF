# -*- coding: utf-8 -*-
"""
Plot time evolution of a field for a given item id
Write a file containing the time evolution of the largest discontinuity
"""

import argparse
import pathlib
import os
import numpy as np
import matplotlib.pyplot as plt
from xfv.src.output_manager.outputdatabaseexploit import OutputDatabaseExploit


def run():
    """
    Run the postprocessing program
    """
    # ----------------------------------------------------------
    # Read instructions
    # ----------------------------------------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("case", help="the path to the output repository")
    parser.add_argument("cut", help="the value of maximal discontinuity opening")
    parser.add_argument("--output_filename", default="all_fields.hdf5",
                        help="the name of the output hdf5 band (default = all_fields.hdf5)")
    parser.add_argument("--save", action="store_true",
                        help="Save the biggest discontinuity opening evolution")
    parser.add_argument("--plot", action="store_true",
                        help="Plot the figure")
    args = parser.parse_args()

    # ----------------------------------------------------------
    # Prepare figure
    # ----------------------------------------------------------
    plt.figure(1)
    plt.xlabel("Time [us]")
    plt.ylabel("DiscontinuityOpening [mm]")

    # ----------------------------------------------------------
    # Plot field evolution for each case
    # ----------------------------------------------------------
    path_to_db = pathlib.Path.cwd().joinpath("..", "tests", args.case, args.output_filename)
    my_hd = OutputDatabaseExploit(path_to_db)

    dict_opening = {}
    ouverture_min = 1.
    ouverture_max = -1.
    the_max_key = -1

    for t in my_hd.saved_times:
        found, disc_opening_at_time = my_hd.extract_discontinuity_opening(t)

        # Fill cell by cell the evolution of discontinuities opening
        if found:
            for i in range(0, len(disc_opening_at_time)):
                index: float = disc_opening_at_time[i, 0]
                if str(index) in dict_opening:  # disc already exists
                    dict_opening[str(index)].append([t, disc_opening_at_time[i, 1]])
                else:
                    dict_opening[str(index)] = [[t, disc_opening_at_time[i, 1]]]

    # Extract relevant info (min and max opening)
    for key in dict_opening:
        array = np.array(dict_opening[key])
        ouverture_min_for_disc = np.min(array[:, 1])
        ouverture_max_for_disc = np.max(array[:, 1])
        ouverture_min = min(ouverture_min, ouverture_min_for_disc)
        if ouverture_max < ouverture_max_for_disc:
            ouverture_max = ouverture_max_for_disc
            the_max_key = key
    print("Ouverture min = " + str(ouverture_min))
    print("Ouverture max =" + str(ouverture_max) + " => cell " + str(the_max_key))

    if args.save:
        np.savetxt("//home/goreckim/Documents/max_opening_{:}.dat".format(
            os.path.basename(args.case)), dict_opening[the_max_key])

    # Plot :
    if args.plot:
        for key in dict_opening:
            array = np.array(dict_opening[key])
            time = array[:, 0] * 1.e6
            mask = (array[:, 1] <= float(args.cut))
            if float(args.cut) > 1e-4:
                plt.plot(time[mask], array[mask, 1] * 1.e3, '.-')
            else:
                plt.plot(time[mask], array[mask, 1] * 1.e6, '.-')
                plt.ylabel("DiscontinuityOpening [um]")

            # limite de critical separation
            if float(args.cut) > 1e-4:
                critical = np.array([[time[mask][0], 1.e-1], [time[mask][-1], 1.e-1]])  # 1.e-4 * 1e3
                plt.plot(critical[:, 0], critical[:, 1], '--', color="black", linewidth=0.1)
        plt.show()


if __name__ == "__main__":
    run()
