#!/usr/bin/env python2.7
# -*- coding: iso-8859-15 -*-
"""
Compare the profil evolution of nodes and cell for xfem simulation with reference
Parameters :
    USER INPUT :
    analysis : type of analysis to be performed : profil ou diff_space

    field : field to be plotted : "VelocityField", "PressureField"

   case_list : list of cases to be considered.
            case : tuple containing : case name, simulation,
                    directory, label and colors pour les graphics
            entrer cette liste entre guillemet et sans espace. séparateur = ','.
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib2tikz import save as tikz_save

import xfv.src.figure_manager.time_figure_tools as fig_tools_path
from xfv.src.utilities.case_definition import CaseManager
from xfv.src.output_figure.profile_tools import get_error_value, plot_field_with_rectangle_shapes, read_hdf5_file, \
    initialize_profile_figure
from xfv.src.output_manager.outputdatabaseexploit import OutputDatabaseExploit

# -----------------------------------------
# Read user instructions
# -----------------------------------------
msg = "Petit programme tout mignon pour tracer les différentes sources de dissipation d'énergie \n"
msg += "Ce script attend comme  arguments  : \n"
msg += "- le type de posttraitement attendu (profile ou diff_space)\n"
msg += "- les champs à tracer (ex: PressureField, VelocityField, DensityField)\n"
msg += "- le temps à regarder (en microsecondes)"
msg += "- une liste des cas à traiter (séparés par une virgule, sans espace, None par défaut) \n"
msg += "- -h ou --help pour afficher l'aide\n"
msg += "Pour une analyse error_profile, le premier cas donné est considéré comme la référence\n"

if sys.argv[1] in ["-h", "--help"]:
    print(msg)
    exit(0)

if len(sys.argv) != 5:
    print(msg)
    exit(0)

# check type of analysis value
analysis = sys.argv[1]
if analysis not in ['profile', 'error_profile']:
    raise ValueError("Type d'analysis {:s} inconnu! Le type doit être parmi profile / diff_space".format(analysis))

# Field to be analyzed
time_field = sys.argv[2].split(',')
time_field = [getattr(fig_tools_path, t_f) for t_f in time_field]
field_list = [t_f.title for t_f in time_field]  # transformation en nom de champ compris par le lecteur de bande hdf5

# Time to plot the profile at
try:
    time_list = sys.argv[3].split(',')
    time_list = [float(t) * 1.e-06 for t in time_list]
except TypeError:
    raise ValueError("Times for profile analysis must be floats")

# Cases of interest
case_list = sys.argv[4].split(',')
if analysis == "error_profile":
    # le premier cas donné est considéré comme la référence
    case_ref = case_list[0]
    case_list = case_list[1:]
    print("Reference case is {:}".format(case_ref.case_name))


# -----------------------------------------
# Initialize treatment and legend
# -----------------------------------------
plt.clf()
plt.close()

draw_discretization_rectangle = False  # représentation du profil avec les histogrammes
plot_all_times_on_same_figure = True  # superpose les courbes pour plusieurs temps
adimensionne_erreur = False  # adimensionne l'erreur commise par le max du choc

# -----------------------------------------
# Run user instruction (analysis)
# -----------------------------------------
for my_case in case_list:
    id_fig = 0
    case = CaseManager().find_case(my_case)
    path_to_db = os.path.join(case.directory_name, "all_fields.hdf5")

    for time in time_list:
        print("------- Infos time = " + str(time) + "------")
        for field in field_list:
            id_fig += 1

            # Creation de la figure :
            if plot_all_times_on_same_figure:
                id_fig = divmod(id_fig, len(field_list))[1]+1
                # title = "{:} field".format(field)
                title = "{:} field at time {:} $\mu s$".format(field, time * 1.e+06)
            else:
                title = "{:} field at time {:} $\mu s$".format(field, time * 1.e+06)
            fig = initialize_profile_figure(id_fig, title)
            plt.ylabel(time_field[divmod(id_fig, len(field_list))[1]-1].label, fontsize=18)
            # Get field index for label if several plots

            # Read hdf5 band :
            my_hd = OutputDatabaseExploit(path_to_db)
            coord, field_value = read_hdf5_file(my_hd, field, time)

            # -----------------------------------------------------------------------------------------
            # Tracé du profil de champ au temps t
            # -----------------------------------------------------------------------------------------
            if analysis == 'profile':
                ax = fig.add_subplot(1, 1, 1)
                # ax.plot(coord * 1.e+03, field_value,
                #         label=case.label, marker=case.marker, linestyle=case.linestyle,
                #         color=case.color)
                ax.plot(coord * 1.e+03, field_value * 1e-9,
                        label=time, marker=case.marker, linestyle=case.linestyle)

                # Plot fields with rectangles showing discretization
                if draw_discretization_rectangle:
                    cell_size = my_hd.get_cells_true_size_at_time(time)
                    plot_field_with_rectangle_shapes(coord, cell_size, field_value, ax)

                ax.legend(loc='best')

            # -----------------------------------------------------------------------------------------
            # Plot the error profile (compared to the chosen reference) at time t
            #  -----------------------------------------------------------------------------------------
            if analysis == 'error_profile':
                # -------------------------------------
                # Read data
                # -------------------------------------
                # -------- Reference field ---------------------------
                path_hd_ref = os.path.join("./0_{:}/".format(case_ref.simulation.upper()),
                                           case_ref.directory_name, "all_fields.hdf5")
                my_hd_ref = OutputDatabaseExploit(path_hd_ref)
                coord_ref, field_value_ref = read_hdf5_file(my_hd_ref, field, time)
                # Initial values
                coord_init, field_value_init = read_hdf5_file(my_hd_ref, field, 0.)

                # Delete last node data in reference because this node does not correspond to any true node in xfv bar
                # (mandatory to get equivalent array size between xfv_left and reference
                if field in ['NodeVelocity', 'NodeCoordinates']:
                    coord_ref = coord_ref[:-1]
                    field_value_ref = field_value_ref[:-1]
                    coord_init = coord_init[:-1]
                    field_value_init = field_value_init[:-1]

                # -------- XVOF field -------------------------------
                enriched_bool = my_hd.extract_field_at_time("CellStatus", time)
                # Extract data of cells at the "left" of discontinuity
                ind_gauche = int(np.where(enriched_bool)[0])
                coord_gauche = coord[:ind_gauche + 1]
                field_gauche = field_value[:ind_gauche + 1]

                # -------------------------------------
                # Compute error
                # -------------------------------------
                coord_interp, error_interp = get_error_value(coord_ref, field_value_ref, coord_gauche, field_gauche,
                                                             error_calcul='absolute')
                erreur_interp_adim = error_interp
                plt.ylabel("Erreur ")
                # Adimensionnement de l'erreur :
                if adimensionne_erreur:
                    max_ref = max(abs(field_value_ref - field_value_init))
                    # correspond à l'amplitude maximale du choc dans le profil
                    print("Adimensionnement de l'erreur par la valeur constante {:} : {:}".format(field, max_ref))
                    erreur_interp_adim = error_interp / max_ref
                    plt.ylabel("Erreur adimensionnee avec le max du choc")

                # Tracé des résultats:
                ax = fig.add_subplot(1, 1, 1)
                ax.plot(coord_interp * 1.e+03, erreur_interp_adim, marker='+',
                        label="Case {:} time {:} $\mu s$".format(case.label, time * 1.e+06))
                ax.legend(loc='best')

tikz_save(os.path.join("//home/goreckim/Documents/These/Matplotlib2Tikz/profil_pression.tex"))
plt.show()
