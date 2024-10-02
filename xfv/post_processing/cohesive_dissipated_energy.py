# -*- coding: utf-8 -*-
"""
Plot the dissipated energy by the cohesive law eventually with experimental data
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
    parser.add_argument("--write_data", action="store_true",
                        help="Write a file with time and dissipated energy")
    args = parser.parse_args()

    if args.case is None:
        raise ValueError("At least one case is needed. Use -case to specify cases to plot")

    # ----------------------------------------------------------
    # Prepare figure
    # ----------------------------------------------------------
    plt.figure(1)
    plt.title("Evolution of the dissipated energy by the cohesive model", fontweight='bold', fontsize=18)
    plt.xlabel("Time [mus]", fontsize=16)
    plt.ylabel("Dissipated energy [J ou J/kg selon les cas]", fontsize=16)


    # ----------------------------------------------------------
    # Read dissipated energy for each cell
    # ----------------------------------------------------------
    for case in args.case[0]:
        dissipated_energy_dict = {}
        path_to_db = pathlib.Path.cwd().joinpath(case, args.output_filename)
        # Read database :
        hd_band = OutputDatabaseExploit(path_to_db)
        for time in hd_band.saved_times:
            cell_status = hd_band.extract_field_at_time("CellStatus", time)[:]
            enriched_cells = np.where(cell_status)[0]
            if len(enriched_cells) > 0:
                dissipated_energy = hd_band.extract_field_at_time("AdditionalDissipatedEnergy", time)[:]
                #print('time=',time)

                for i in range(len(dissipated_energy[:, 0])):
                    cell_id = dissipated_energy[i, 0]
                    dis_energy = dissipated_energy[i, 1]
                    try:
                        dissipated_energy_dict[cell_id].append([time, dis_energy])
                    except KeyError:  # 1ere fois que la maille est enrichie => init list opening
                        dissipated_energy_dict[cell_id] = [[time, dis_energy]]
                    #print(dissipated_energy_dict[cell_id][0,1])

        # ----------------------------------------------------------
        # Permet d'aficher les energies dissipées pour toutes les cellules
        # ----------------------------------------------------------         
        if args.item_ids[0][0] == 1000:
            args.item_ids[0] = []
            args.item_ids[0] = dissipated_energy[:,0]
        

        #Transformation en array
        for key in dissipated_energy_dict:
            dissipated_energy_dict[key] = np.array(dissipated_energy_dict[key])

        # import ipdb; ipdb.set_trace()

        for item_id in args.item_ids[0]:
            # Read database :
            j = 0.
            for ane in dissipated_energy[:, 0]:
                if ane == item_id:
                    j += 1.
            if j < 0.99:
                print(dissipated_energy[:, 0])
                exit(f'item_ids selectionné ne fait pas parti de la liste des cellules enrichies. Choississez item_ids (ie. le numéro de cellules) dans la liste ci-dessus ou tapez 1000 pour selectionner toutes les cellules')
            # Plot field :
            plt.plot(dissipated_energy_dict[item_id][:,0],dissipated_energy_dict[item_id][:,1],label= 'cell n°'+str(item_id))
            if args.write_data:
                data_path = f"Field_evolution_cohesive_dissipated_energy_at_cell_{item_id}.dat"
                with open(data_path, "w") as file_object:
                    for x_data, y_data in zip(dissipated_energy_dict[item_id][:,0], dissipated_energy_dict[item_id][:,1]):
                        file_object.write("{:20.18g}\t{:20.18g}\n".format(x_data, y_data))
                print("Data written in {:s}".format(data_path))

if __name__ == "__main__":
    run()
    # ----------------------------------------------------------
    # Show figure
    # ----------------------------------------------------------
    plt.legend(loc="right")
    plt.show()
