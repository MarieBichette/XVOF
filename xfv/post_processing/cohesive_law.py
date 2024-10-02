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
    parser.add_argument("-item_ids", type=int, action='append', nargs='+', help="the id of the item to look at")
    parser.add_argument("--output_filename", default="all_fields.hdf5",
                        help="the name of the output hdf5 band (default = all_fields.hdf5)")
    parser.add_argument("--write_data", action="store_true", help="To write data in an output file")
    args = parser.parse_args()

    if args.case is None:
        raise ValueError("At least one case is needed. Use -case to specify cases to plot")

    # ----------------------------------------------------------
    # Prepare figure
    # ----------------------------------------------------------
    plt.figure(1)
    plt.title("Evolution of the cohesive law", fontweight='bold', fontsize=18)
    plt.xlabel("discontinuity opening [m]", fontsize=16)
    plt.ylabel("cohesive force []", fontsize=16)


    # ----------------------------------------------------------
    # Read discontinuity opening for each cell
    # ----------------------------------------------------------
    for case in args.case[0]:
        opening_dict = {}
        force_dict = {}
        path_to_db = pathlib.Path.cwd().joinpath(case, args.output_filename)
        # Read database :
        hd_band = OutputDatabaseExploit(path_to_db)
        for time in hd_band.saved_times:
            cell_status = hd_band.extract_field_at_time("CellStatus", time)[:]
            enriched_cells = np.where(cell_status)[0]
            if len(enriched_cells) > 0:
                opening = hd_band.extract_field_at_time("AdditionalDiscontinuityOpening", time)[:]
                cohesive_force = hd_band.extract_field_at_time("AdditionalCohesiveForce", time)[:]
                for i in range(len(opening[:, 0])):
                    cell_id = opening[i, 0]
                    op = opening[i, 1]
                    force = cohesive_force[i, 1]
                    try:
                        opening_dict[cell_id].append([time, op])
                        force_dict[cell_id].append([time, force])
                    except KeyError:  # 1ere fois que la maille est enrichie => init list opening
                        opening_dict[cell_id] = [[time, op]]
                        force_dict[cell_id] = [[time, force]]
        
    # Permet de selectionner les item_ids pour toutes les cellules enrichies
    if args.item_ids[0][0] == 1000:
        args.item_ids[0] = []
        args.item_ids[0] = opening[:,0]
        print(args.item_ids[0])

    #Transformation en array
    for key in opening_dict:
        opening_dict[key] = np.array(opening_dict[key])
    for key in force_dict:
        force_dict[key] = np.array(force_dict[key])
        # import ipdb; ipdb.set_trace()

    for item_id in args.item_ids[0]:
        # Read database :
        j = 0.
        for ane in opening[:, 0]:
            if ane == item_id:
                j += 1.
                if j < 0.99:
                    print(opening[:, 0])
                    exit(f'Choissisez item_ids dans la liste ci-dessus')
            # Plot field :
        plt.plot(opening_dict[item_id][:,1],force_dict[item_id][:,1], '+', label= 'cell n°'+str(item_id))

        if args.write_data:
            data_path = f"Field_evolution_cohesive_law_{item_id}.dat"
            with open(data_path, "w") as file_object:
                for x_data, y_data in zip(opening_dict[item_id][:,1], force_dict[item_id][:,1]):
                    file_object.write("{:20.18g}\t{:20.18g}\n".format(x_data, y_data))
            print("Data written in {:s}".format(data_path))

if __name__ == "__main__":
    run()
    # ----------------------------------------------------------
    # Show figure
    # ----------------------------------------------------------
    plt.legend(loc="best")
    plt.show()
