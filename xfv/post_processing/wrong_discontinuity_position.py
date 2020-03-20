#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
Script to characterize the error in thermodynamics when the discontinuity is artificially placed in the cracked element.
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

import xvof.src.figure_manager.time_figure_tools as fig_tools_path
from xvof.src.utilities.case_definition import CaseManager
from xvof.src.output_figure.hdf5_posttreatment_tools import read_database_in_array
from xvof.src.output_manager.outputdatabaseexploit import OutputDatabaseExploit


print("Ce script estamené à disparaitre. Il faut utiliser show_error_with_epsilon")
os._exit(0)

try:
    field = sys.argv[1]
    # recognition of the string field
    field = getattr(fig_tools_path, field)
except IndexError:
    raise ValueError("""Ce script attend comme argument le nom du champ à tracer""")


example_case = CaseManager().find_case("eps_002_hansbo_lump_somme")


# Récupération des données quand on impose un position de la discontinuité
imposed_epsilon = np.array([0.1])
my_hd = OutputDatabaseExploit(os.path.join(example_case.directory_name, "all_fields.hdf5"))
size = len(my_hd.saved_times)
imposed_field_cell_500 = np.zeros([size, len(imposed_epsilon)])
case_number = 0


for imposed_eps in imposed_epsilon:
    if imposed_eps != 0.1:
        imposed_case = CaseManager().find_case("eps_00{:}_hansbo_lump_somme".format(int(imposed_eps*100)))
    else:
        imposed_case = CaseManager().find_case("eps_010_hansbo_lump_somme")
    print "Imposed epsilon case : "
    res = read_database_in_array(os.path.join(imposed_case.directory_name, "all_fields.hdf5"), 500, field.title)[:, 1]
    imposed_field_cell_500[:, case_number] = res
    case_number += 1

# Récupération des données quand on utlise la vraie position de la discontinuité
nbr_of_real_position = 10
real_field_in_cell_500 = np.zeros([size, nbr_of_real_position])
real_epsilon = np.zeros([nbr_of_real_position])
# for i in range(nbr_of_real_position):
for i in [1,3,5,7,9]:
    epsilon = (i +1) / 100.
    print "Real case " + str(epsilon)
    real_epsilon[i] = epsilon
    if epsilon != 0.1:
        current_case = CaseManager().find_case("eps_00{:}_hansbo_lump_somme".format(i+1))
    else:
        current_case = CaseManager().find_case("eps_010_hansbo_lump_somme")
    print current_case.directory_name
    res = read_database_in_array(os.path.join(current_case.directory_name, "all_fields.hdf5"), 500, field.title)[:, 1]
    real_field_in_cell_500[:, i] = res

fig = plt.figure(99)
fig.patch.set_color("white")
fig.suptitle("Error in {:} with arbitrary $\epsilon$ (%)". format(field.title), fontsize=20, fontweight='bold')
plt.xlabel("Discontinuity relative position [-]", fontsize=20)
plt.ylabel("Error in pressure [%]", fontsize=20)


erreur_epsilon = np.zeros([nbr_of_real_position, len(imposed_epsilon)])
my_array = np.zeros([my_hd.nb_saved_times, len(real_epsilon), len(imposed_epsilon)])

for i_imp in range(len(imposed_epsilon)):
    for j_real in range(len(real_epsilon)):
        my_array[:, j_real, i_imp] = (imposed_field_cell_500[:, i_imp] - real_field_in_cell_500[:, j_real]) / \
                                        max(np.abs(imposed_field_cell_500[:, i_imp])) * 100

erreur_field = np.max(my_array, 0)

for i in range(len(imposed_epsilon)):
    plt.plot(real_epsilon, erreur_field[:, i], '*-', label="arbitrary $\epsilon$ = {:}".format(imposed_epsilon[i]))


plt.legend(loc='best', fontsize=15)
plt.show()
