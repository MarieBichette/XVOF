#!/usr/bin/env python2.7
# -*- coding:  iso-8859-15 -*-
"""
Figure divergence de la vitesse en fonction du temps
Paramètres d'entrée :
- cas à analyser
- type d'analyse : divergence ou vitesse du bord de la fissure reconstruite
"""

import os
import sys
from matplotlib2tikz import save as tikz_save
import matplotlib.pyplot as plt
import numpy as np

from xfv.src.utilities.case_definition import CaseManager
from xfv.src.output_manager.outputdatabaseexploit import OutputDatabaseExploit
from xfv.src.utilities.velocity_data_posttraitment import compute_velocity_divergence_at_time, \
    compute_velocity_on_discontinuity_borders_at_time

msg = "Mini programme pour tracer la divergence du champ de vitesse ou le champ de vitesse sur les bords de la disc\n" \
      "Ce script attend comme arguments  : \n"
msg += "- le type de tracé (divergence | bordervelocity) \n"
msg += "- la list des cas à traiter \n"
msg += "- -h ou --help pour afficher l'aide\n"

if sys.argv[1] in ["-h", "--help"]:
    print(msg)
    exit(0)

if len(sys.argv) != 3:
    print(msg)
    exit(0)

project_dir = os.path.split(os.path.dirname(os.path.abspath(os.path.curdir)))[0]
base_dir = os.path.split(project_dir)[0]

case_list = sys.argv[2].split(',')
analysis = sys.argv[1]
if analysis not in ["divergence", "bordervelocity"]:
    raise ValueError (""" Enter analysis type in (divergence | bordervelocity)""")

# Creation de la figure :
fig = plt.figure(1)
fig.patch.set_facecolor("white")
plt.xlabel("Time [$\mu s$]", fontsize=18)

if analysis == "divergence":
    fig.suptitle("Time evolution of velocity field divergence", fontsize=14, fontweight='bold')
    plt.ylabel("div $u$ [$s^{-1}$]", fontsize=18)
else:
    fig.suptitle("Time evolution of boundary velocity of the discontinuity", fontsize=14, fontweight='bold')
    plt.ylabel("Velocity [$ms^{-1}$]", fontsize=18)

# -----------------------------------------
# Run user instruction (analysis)
# -----------------------------------------
for case in case_list:
    print("---------------------------")
    case = CaseManager().find_case(case)
    print(case.case_name)
    path_to_db = os.path.join(case.directory_name, "all_fields.hdf5")
    print(path_to_db)
    my_hd = OutputDatabaseExploit(path_to_db)

    time = np.array(my_hd.saved_times)
    time_index = 0
    field_to_plot = np.zeros(my_hd.nb_saved_times)

    for t in my_hd.saved_times:
        if analysis == "divergence":
            # import ipdb ; ipdb.set_trace()
            field_to_plot[time_index] = compute_velocity_divergence_at_time(my_hd, t, 500)
        else:
            field_to_plot[time_index] = compute_velocity_on_discontinuity_borders_at_time(my_hd, t, [501])
        time_index += 1

    # -----------------------------------------------------------------------------------------
    # Tracé du profil de champ au temps t
    # -----------------------------------------------------------------------------------------
    ax = fig.add_subplot(111)

    ax.plot(time * 1.e+06, field_to_plot,
            marker=case.marker, markersize=5, linestyle=case.linestyle, color=case.color, label=case.label)

    if analysis=="divergence":
        ax.set_xlim([2.5, 5.])
        ax.set_ylim([-15000, 10000])

plt.legend(loc='best')

tikz_save(os.path.join(base_dir, "Documents/These/Matplotlib2Tikz/velocity_divergence.tex"))

plt.show()
