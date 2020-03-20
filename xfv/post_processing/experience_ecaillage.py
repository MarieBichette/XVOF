#!/usr/bin/env python2.7
# -*- coding: iso-8859-15 -*-
""" 
Compare velocities of free surfaces from simulation with experimental data
USER INPUT :
- The case name to be compared with experiment
(field to be compare is VelocityField for the moment)
"""

import os
import sys
import matplotlib.pyplot as plt
import numpy as np

from xfv.src.utilities.case_definition import CaseManager
from xfv.src.figure_manager.time_figure_tools import VelocityField
from xfv.src.output_figure.hdf5_posttreatment_tools import read_database_in_array
from xfv.src.utilities.experimental_data_exploit import ExperimentalData

# -----------------------------------------
# Initialize treatment and legend
# -----------------------------------------
plt.clf()
plt.close()


project_dir = os.path.split(os.path.dirname(os.path.abspath(os.path.curdir)))[0]
base_dir = os.path.split(project_dir)[0]

english = False
extension = '.hdf5'

recalage_temps = True
id_item = -1  # dernier noeud (surface libre du bord opposé)

experimental_file = "TA6.TXT"
experimental_dir = os.path.join(base_dir, "Documents/These/Experimental_data/Ta/")

exp_data = os.path.join(experimental_dir, experimental_file)

hdf5_band = "all_fields.hdf5"

fig = plt.figure(VelocityField.colonne_history)
fig.patch.set_facecolor("white")
if english:
    fig.suptitle("Evolution of the free surface velocity", fontsize=20, fontweight='bold')
    plt.xlabel("Time [$\mu$s]")
    plt.ylabel("Velocity of free surface [m/s]")
else:
    fig.suptitle("Vitesse de surface libre", fontsize=20, fontweight='bold')
    # plt.xlabel("Temps [$\mus$]", fontsize=18)
    plt.ylabel("Vitesse de la surface libre [m/s]", fontsize=18)

# -----------------------------------------
# Read user instructions
# -----------------------------------------
msg = "Petit programme tout mignon pour tracer la vitesse de surface libre simulée et expérimentale \n"
msg += "Ce script attend comme  arguments  : \n"
msg += "- une liste des cas à traiter (séparés par une virgule, sans espace, None par défaut) \n"
msg += "- -h ou --help pour afficher l'aide\n"
msg += "On pose la question si on souhaite afficher les données expérimentales actuelles " \
       "(adresse du fichier dans le code)"

if len(sys.argv) != 2:
    print msg
    exit(0)

if sys.argv[1] in ["-h", "--help"]:
    print msg
    exit(0)

# choix_exp = raw_input("Voulez vous tracer les résultats expérimentaux du fichier {:} ? (y/n) ".format(exp_data))
choix_exp = "y"

case_list = sys.argv[1].split(',')
case_list = [CaseManager().find_case(case) for case in case_list]

# Define and plot experimental results
if choix_exp == "y":
    exp_result = ExperimentalData(exp_data)
    exp_result.compute_spall_stress(rho_0=16960., c_0=3460., us1=274, us2=109)
    # exp_result.compute_spall_stress(rho_0=8930., c_0=3940., us1=502, us2=431)
    exp_result.plot_experimental_results(VelocityField.colonne_history)

# -----------------------------------------
# Run user instruction (analysis)
# -----------------------------------------
for case in case_list:
    print "-------------------------------"
    print case.case_name
    # exp_result.plot_hugoniot(case)

    # Tracé des données numériques
    path_to_hdf5_db = os.path.join(case.directory_name, hdf5_band)
    print "Read hdf5 data in file {:}".format(path_to_hdf5_db)

    hdf5_data = read_database_in_array(path_to_hdf5_db, id_item, "NodeVelocity", 1)
    time = hdf5_data[:, 0]
    velocity = hdf5_data[:, 1]
    mask = (velocity > 1)
    time_0 = time[np.where(mask)[0][0]]

    np.savetxt(os.path.join(case.directory_name, "vitesse_surface_libre.dat"), hdf5_data)

    plt.figure(VelocityField.colonne_history)
    case_label = u"{:}".format(case.label)
    plt.plot((time - time_0) * 1.e+06, velocity, color=case.color, marker=case.marker, label=case_label, linestyle=case.linestyle)
    # plt.plot((time[mask] - time_0) * 1.e+06, velocity[mask], color=case.color, marker=case.marker, label=case.label)

plt.legend(loc="best")
plt.show()
