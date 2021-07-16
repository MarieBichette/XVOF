# -*- coding: utf-8 -*-
"""
Plot the free surface velocity eventually with experimental data
"""

import argparse
import pathlib
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple

from xfv.src.output_manager.outputdatabaseexploit import OutputDatabaseExploit

CellInfo = namedtuple("CellInfo", ["ouverture_min", "ouverture_max", "temps_apparition"])

def run():
    """
    Run post processing program
    """
    # ----------------------------------------------------------
    # Read instructions
    # ----------------------------------------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("-case", action='append', nargs='+',
                        help="the path to the output repository")
    parser.add_argument("--output_filename", default="all_fields.hdf5",
                        help="the name of the output hdf5 band (default = all_fields.hdf5)")
    args = parser.parse_args()

    if args.case is None:
        raise ValueError("At least one case is needed. Use -case to specify cases to plot")

    # ----------------------------------------------------------
    # Prepare figure
    # ----------------------------------------------------------
    plt.figure(1)
    plt.title("Evolution of the free surface velocity", fontweight='bold', fontsize=18)
    plt.xlabel("Time [mus]", fontsize=16)
    plt.ylabel("Free surface velocity [m/s]", fontsize=16)


    # ----------------------------------------------------------
    # Read discontinuity opening for each cell
    # ----------------------------------------------------------
    for case in args.case[0]:
        opening_dict = {}
        path_to_db = pathlib.Path.cwd().joinpath(case, args.output_filename)
        # Read database :
        hd_band = OutputDatabaseExploit(path_to_db)
        for time in hd_band.saved_times:
            cell_status = hd_band.extract_field_at_time("CellStatus", time)[:]
            enriched_cells = np.where(cell_status)[0]
            if len(enriched_cells) > 0:
                opening = hd_band.extract_field_at_time("AdditionalDiscontinuityOpening", time)[:]
                for i in range(len(opening[:, 0])):
                    cell_id = opening[i, 0]
                    op = opening[i, 1]
                    try:
                        opening_dict[cell_id].append([time, op])
                    except KeyError:  # 1ere fois que la maille est enrichie => init list opening
                        opening_dict[cell_id] = [[time, op]]

        #Transformation en array
        for key in opening_dict:
            opening_dict[key] = np.array(opening_dict[key])

        # import ipdb; ipdb.set_trace()

        # Ouverture min / max :

        print(case + " : " + str(len(opening_dict)) + " disc créées")


if __name__ == "__main__":
    run()
    # ----------------------------------------------------------
    # Show figure
    # ----------------------------------------------------------
    # plt.legend(loc="best")
    # plt.show()
