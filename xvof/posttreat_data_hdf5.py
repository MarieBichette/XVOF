#!/usr/bin/env python2.7
#  -*- coding: iso-8859-1 -*-
"""
Posttreat the hdf5 band to write the true fields for cells and nodes selected
Call with arg :
item : cell/node
id_item : int (list of int possible , entre guillemet et sans espace)
simulation : ref, xfem
case (optional) : list
"""

import os
import sys
import numpy as np
from xvof.output_manager.outputdatabaseexploit import OutputDatabaseExploit
from xvof.figure_manager.time_figure_manager import TimeFigureTools
import xvof.case_definition as case_path
from collections import namedtuple
# -----------------------------------------
# Read user instructions
# -----------------------------------------
if len(sys.argv) > 5:
    raise ValueError("Ce script attend comme  arguments  : \nle type d'item (node|cell) \n le numéro de l'item "
                     "(de type int) \nle type de simulation (ref|xfem) \net une liste des cas à traiter"
                     " (case_list : liste des cas à traiter (entre guillemet et sans espace) | None par défaut) ")
item = sys.argv[1]
if item not in ['node', 'cell']:
    raise ValueError("Type d'item {:s} inconnu! Le type doit être soit node soit cell".format(item))

id_item_list_str = sys.argv[2].split(',')
try:
    id_item_list = map(int, id_item_list_str)
except:
    raise ValueError("Le numéro d'item doit être un entier ou liste d'entiers")

try:
    case_list = sys.argv[3].split(',')
except ValueError:
    case_list = [case_path.default]
    print "No case specified. Taking default database stored in path 0_XFEM / 0_REF\n" \
          "Enter a list of cases entre guillemet et sans espace"  # à améliorer...
# -----------------------------------------
# Initialize treatment and legend
# -----------------------------------------
if case_list == "all_mass":
    case_list = ["mass_1", "mass_2", "mass_3", "mass_4", "mass_5", "mass_6", "mass_7", "mass_8"]
for my_case in case_list:
    case = getattr(case_path, my_case)
    if case.case_name == "analytical":
        case_list[case_list.index(case)] = "analytical_inverse_complete_matrix"

my_dico_array_to_write = dict()
# -----------------------------------------
# Read and save databases for selected cases
# -----------------------------------------
for my_case in case_list:
    case = getattr(case_path, my_case)
    simulation = case.simulation
    my_figure_tools = TimeFigureTools(case, item, id_item_list[0])
    my_figure_tools.find_field_list()

    # Find database location:
    # import ipdb; ipdb.set_trace()
    directory = "./0_{:s}/{:}".format(simulation.upper(), case.directory_name)
    path_to_hdf5_db = os.path.join(directory, "all_fields.hdf5")

    # Read data base and post treat
    my_hd = OutputDatabaseExploit(path_to_hdf5_db)
    for id_item in id_item_list:
        my_dico_array_to_write[str(id_item)] = np.zeros([my_hd.nb_saved_times, len(my_figure_tools.field_list) + 1])

    indice_temps = 0
    for t in my_hd.saved_times:
        time_registration_to_be_done = True
        for field in my_figure_tools.field_list:
            for id_item in id_item_list:
                if time_registration_to_be_done: # On enregistre qu'une seule fois chaque temps par champs
                    my_dico_array_to_write[str(id_item)][indice_temps, 0] = t
                my_dico_array_to_write[str(id_item)][indice_temps, field.colonne_history] = \
                    my_hd.extract_true_field_at_time(field.title, t)[id_item][1]
            time_registration_to_be_done = False
        indice_temps += 1

    # Save cell or node fields in appropriate .dat file
    for id_item in id_item_list:
        header = "-----Time History for {:} {:} :".format(item, str(id_item)) + os.linesep
        path_to_save_data = os.path.join(directory, "{:s}_history_{:<03d}.dat".format(item, id_item))
        with open(path_to_save_data, 'w') as f:
            if item == 'cell':
                np.savetxt(f, my_dico_array_to_write[str(id_item)],
                           fmt=['%+10.9e', '%+10.9e', '%+10.9e', '%+10.9e', '%+10.9e'], header=header)
            if item == 'node':
                np.savetxt(f, my_dico_array_to_write[str(id_item)], fmt=['%+10.9e', '%+10.9e'], header=header)
            f.close()
print "Done !"
# Trop long --> paralléliser les calculs
