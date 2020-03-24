#!/usr/bin/env python2.7
# -*- coding: utf-8-*-

"""
Script pour tracer les fichiers de données sur le modèle cohésif.
"""

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from xfv.src.output_manager.outputdatabaseexploit import OutputDatabaseExploit
from xfv.src.utilities.case_definition import CaseManager


def prepare_figure(fig_id, title, x_label, y_label):
    fig = plt.figure(fig_id)
    fig.patch.set_facecolor("white")
    plt.title(title, fontsize=20, fontweight='bold')
    plt.xlabel(x_label, fontsize=20)
    plt.ylabel(y_label, fontsize=20)

msg = "Petit programme tout mignon pour tracer les lescohésives sur chaque discontinuité \n"
msg += "Ce script attend comme  arguments  : \n"
msg += "- une liste des cas à traiter (séparés par une virgule, sans espace, None par défaut) \n"
msg += "- -h ou --help pour afficher l'aide\n"

if len(sys.argv) != 2:
    print(msg)
    exit(0)

if sys.argv[1] in ["-h", "--help"]:
    print(msg)
    exit(0)

case_list = sys.argv[1].split(',')
case_list = [CaseManager().find_case(case) for case in case_list]

color_list = ['blue','red', 'orange', 'green', 'purple', 'gray', 'skyblue', 'lightgreen', 'darkgreen'] * 100
# prepare_figure(1, "Cohesive laws", "Opening [$\mu m$]", "Cohesive force [$Pa$]")
prepare_figure(1, "Comportement des zones cohesives", "Ouverture $\delta$ [$\mu m$]", "Contrainte cohesive $\sigma$ [$Pa$]")
# prepare_figure(2, "Force vs time", "Temps [s]", "Cohesive force [Pa]")
# prepare_figure(3, "Opening vs time", "Temps [s]", "Discontinuity opening [mm]")


for case in case_list:
    output_db = OutputDatabaseExploit(os.path.join(case.directory_name, "all_fields.hdf5"))
    # nombre max de  discontinuités :
    disc_collection_final = np.where(output_db.extract_field_at_time("CellStatus", output_db.saved_times[-1]))[0]
    nb_disc_total = len(disc_collection_final)

    print("Nombres de  discontinuités créées : {:}".format(nb_disc_total))

    # gros tableau pour les résultats
    opening = np.zeros([output_db.nb_saved_times, nb_disc_total])
    force = np.zeros([output_db.nb_saved_times, nb_disc_total])
    temps = np.zeros([output_db.nb_saved_times, nb_disc_total])
    mask_existance = np.zeros([output_db.nb_saved_times, nb_disc_total], dtype=bool)
    indice_cell_rompue = np.zeros([output_db.nb_saved_times, nb_disc_total])

    # Données pour tracer les axes sur les courbes après
    min_force = 10.e+10
    max_force = 0.
    min_opening = 10.e+10
    max_opening = 0.

    for i_temps in range(output_db.nb_saved_times):
        t = output_db.saved_times[i_temps]
        temps[i_temps, :] = np.ones(nb_disc_total) * t

        for i_disc in range(nb_disc_total):
            exist, op, f = output_db.extract_fields_for_cohesive_zone_model(disc_collection_final[i_disc], t)

            opening[i_temps, i_disc] = op * 1.e+03 # multiplicateur pour tout passer en micro mètre
            force[i_temps, i_disc] = f
            mask_existance[i_temps, i_disc] = exist
            indice_cell_rompue[i_temps, i_disc] = disc_collection_final[i_disc]

        # sauvegarde des données utiles pour le tracé des axes sur les courbes :
        min_force = min(min_force, min(force[i_temps, :]))
        max_force = max(max_force, max(force[i_temps, :]))
        min_opening = min(min_opening, min(opening[i_temps, :]))
        max_opening = max(max_opening, max(opening[i_temps, :]))

    # Tracé des résultats par discontinuité :
    # for disc_j in xrange(nb_disc_total):
    for disc_j in [6,0]:

        # Filtrage des données : on s'intéresse uniquement à la discontinuité une fois qu'elle a été créée
        opening_disc_j = opening[:, disc_j][mask_existance[:, disc_j]]
        force_disc_j = force[:, disc_j][mask_existance[:, disc_j]]
        temps_disc_j = temps[:, disc_j][mask_existance[:, disc_j]]
        indice_cell_rompue_j = int(indice_cell_rompue[:, disc_j][mask_existance[:, disc_j]][0])

        # Trace de la loi cohésive force en fonction de ouverture
        plt.figure(1)
        # plt.plot(opening_disc_j, force_disc_j, '.-',
        #          label="{:} - disc on cell {:}".format(case.label, indice_cell_rompue_j),
        #          color=color_list[disc_j])
        plt.plot(opening_disc_j, force_disc_j, '.-',color=color_list[disc_j])

        # identification du dernier point sur la courbe
        plt.plot(opening_disc_j[-1], force_disc_j[-1], 'd', color=color_list[disc_j])
        plt.plot([0, 0], [min_force, max_force], '--', color="black")
        plt.plot([min_opening, max_opening], [0, 0], '--', color="black")
        # plt.legend(loc='best')

        # Tracé de force en fonction du temps
        # plt.figure(2)
        # plt.plot(temps_disc_j, force_disc_j, '.-',
        #          label="{:} - disc on cell {:}".format(case.label, indice_cell_rompue_j),
        #          color=color_list[disc_j])
        # plt.plot([0, 0], [min_force, max_force], '--', color="black")
        # plt.plot([0, output_db.saved_times[-1]], [0, 0], '--', color="black")
        # plt.legend(loc='best')

        # Tracé de ouveture en fonction du temps
        # plt.figure(3)
        # plt.plot(temps_disc_j, opening_disc_j, '.-',
        #          label="{:} - disc on cell {:}".format(case.label, indice_cell_rompue_j),
        #          color=color_list[disc_j])
        # plt.plot([0, 0], [min_opening, max_opening], '--', color="black")
        # plt.plot([0, output_db.saved_times[-1]], [0, 0], '--', color="black")
        # plt.legend(loc='best')
plt.show()

