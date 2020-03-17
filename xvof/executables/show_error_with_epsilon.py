#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
Script to plot the difference between a reference and enriched case versus the relative discontinuity position
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from xvof.src.utilities.case_definition import CaseManager
from xvof.src.figure_manager.time_figure_manager import TimeFigureManager
from xvof.src.figure_manager.time_figure_tools import TimeFigureTools, DensityField


def translate_epsilon_to_case(case_name_base, epsilon):
    if epsilon < 0.1:
        case = "eps_00{:}".format(int(epsilon * 100))
    else:
        case = "eps_0{:}".format(int(epsilon * 100))
    case += case_name_base
    return CaseManager().find_case(case)


# -----------------------------------------
# Read user instructions
# -----------------------------------------
msg = "Petit programme tout mignon pour tracer les écarts entre ref et enrichissement pour la comparaison hydro \n"
msg += "Ce script attend comme  arguments  : \n"
msg += "- une liste des epsilon à tracer (float séparés par des -)\n"
msg += "Si les données de références ont des epsilons différents : une liste des epsilon pour la référence " \
       "(float séparés par des -) \n"
msg += "On demandera une option pour savoir comment calculer l'erreur (absolue, adimensionnée, relative) \n"
msg += "- -h ou --help pour afficher l'aide"

if len(sys.argv) not in [2, 3]:
    print msg
    exit(0)

##########################
# todo : à coder pour remplacer wrong disc. position
if len(sys.argv) == 3:
    print("pas encore codé")
    # récuéprer la liste des epsilon et init de la liste des ref avec ça
    os._exit(0)
###########################

if sys.argv[1] in ["-h", "--help"]:
    print msg
    exit(0)

try:
    eps_list = [float(eps) for eps in sys.argv[1].split('-')]
except TypeError:
    print msg
    raise

# On demande quel calcul d'erreur faire
option = None
while option not in range(3):
    option = input("Quelle méthode de calcul pour l'erreur ? (0 = absolue, 1 = adimensionnée, 2 = relative) ")

# -----------------------------------------
# Construction de la liste des cas
# ------------------------------------------
case_hansbo_base = "_hansbo"
case_lump_eps_base = "_hansbo_lump_eps"
case_lump_sum_base = "_hansbo_lump_somme"
case_ref_base = "_ref_wilkins"

nb_case_epsilon = len(eps_list)
# Tableau d'erreur max :
error_max = np.zeros([nb_case_epsilon, 3])
# Array hansbo | lump_eps | lump_sum

field = DensityField
print "-> Selected field is : " + field.label

hansbo_cases = [translate_epsilon_to_case(case_hansbo_base, epsilon) for epsilon in eps_list]
hansbo_lump_eps_cases = [translate_epsilon_to_case(case_lump_eps_base, epsilon) for epsilon in eps_list]
hansbo_lump_sum_cases = [translate_epsilon_to_case(case_lump_sum_base, epsilon) for epsilon in eps_list]
enriched_cases = np.array([hansbo_cases, hansbo_lump_eps_cases, hansbo_lump_sum_cases])

ref_cases = [translate_epsilon_to_case(case_ref_base, epsilon) for epsilon in eps_list]
# import ipdb ; ipdb.set_trace()

simulations = ["Hansbo consistant", "Hansbo lump Menouillard", "Hansbo lump somme"]
answer_array = np.zeros([len(simulations)], dtype=bool)

# construction d'un vecteur erreur pour chaque simulation
sim_id = 0
while sim_id < len(simulations):
    answer = raw_input("Voulez vous étudier le cas : {:} ? [y/n]".format(simulations[sim_id]))
    answer_array[sim_id] = (answer == "y")
    if answer == "y":
        for i in range(nb_case_epsilon):
            case_ref = ref_cases[i]
            case_enr = enriched_cases[sim_id, i]

            # datafile = os.path.join(case_enr.directory_name, "all_fields.hdf5")
            # import ipdb; ipdb.set_trace()

            fig_manager = TimeFigureManager(case_enr, "cell", 500, "all_fields.hdf5", case_ref=case_ref)
            error_at_all_t = fig_manager.compute_time_error_single(field, compute_error_index=option)
            error_max[i, sim_id] = max(error_at_all_t)
    sim_id += 1

# tracé des résultats :
fig = plt.figure(1)
fig.patch.set_color("white")
plt.title("Erreur sur la densite pour differents choix de matrice de masse")
plt.xlabel("Epsilon")
plt.ylabel("Erreur maxi sur la densite")
sim_id = 0
while sim_id < len(simulations):
    plt.plot(eps_list, error_max[:, sim_id], label=simulations[sim_id])
    sim_id += 1

plt.legend()
plt.show()
