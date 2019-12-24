#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
A script writing a march diagram after exploitation of the output database
"""
from __future__ import print_function
import os
import sys
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from xvof.utilities.case_definition import CaseManager
from xvof.output_manager.outputdatabaseexploit import OutputDatabaseExploit
from xvof.executables.space_time_diagram_utilities import SpaceTimeDiagramUtilities

# ----------------------------------------------------------------
# Quelques paramètres
# ----------------------------------------------------------------
plot = True
english = False
matplotlib.rcParams['text.usetex'] = False
matplotlib.rcParams['text.latex.unicode'] = False

# ----------------------------------------------------------------
# Définition des labels et légendes
# ----------------------------------------------------------------
if english:
    my_xlabel = "Position [mm]"
    my_ylabel = "Time [$\mu$s]"
    my_base_title = ["Space - Time diagram in Pressure", "Space - Time diagram "]
    my_colorbarlabel = ["Pressure", "Pressure gradient"]
else:
    my_xlabel = "Position [mm]"
    my_ylabel = "Temps [$\mu$s]"
    my_base_title = ["Diagramme de marche en Pression", "Diagramme de marche"]
    my_colorbarlabel = ["Pression", "Gradient de pression"]

# ----------------------------------------------------------------
# Création d'une colormap jolie
# ----------------------------------------------------------------
cdict = {'red': ((0., 0., 0.),
                  (0.5, 0.75, 0.75),
                  (1., 1., 1.)),
          'green': ((0., 0., 0.),
                  (0.5, 0.75, 0.75),
                  (1., 0., 0.)),
            'blue': ((0., 0.5, 0.5),
                  (0.5, 0.75, 0.75),
                  (1., 0., 0.))}
my_cmap = LinearSegmentedColormap('custom_cmap', cdict)
plt.register_cmap(cmap=my_cmap)

# ----------------------------------------------------------------
# Lecture des instructions utilisateur
# ----------------------------------------------------------------
msg = "Mini programme pour tracer un diagramme de marche\n" \
      "Ce script attend comme arguments  : \n"
msg += "- le type de simulation (classic|schlieren|coupe|pick) \n"
msg += "- le type de champ à tracer (ex : Pressure)\n"
msg += "- les cas à traiter (séparés par une virgule, sans espace)\n"
msg += "- -h ou --help pour afficher l'aide\n"

# Affichage de l'aide
if len(sys.argv) != 4:
    print(msg)
    exit(0)

if sys.argv[1] in ["-h", "--help"]:
    print(msg)
    exit(0)

# Type d'analyse
analysis = sys.argv[1]
analysis_int = 0
if analysis not in ["schlieren", "classic", "coupe", "pick"]:
    print("Type d'analyse inconnu. Choisir entre classic et schlieren et coupe et pick")
    exit(0)
if analysis == 'schlieren':
    analysis_int = 1

# Demande à l'utilisateur quelle coupe il veut tracer
if analysis == "coupe":
    coupe_type = -1
    while coupe_type not in [0, 1]:
        coupe_type = int(raw_input("Quel type de coupe fait-on ? (0 = constant time, 1 = constant item) << "))

# Cas demandé
case_list = sys.argv[3].split(',')
case_list = [CaseManager().find_case(my_case) for my_case in case_list]

# Champ demandé (la vérification se fait après avoir lu les champs qui existent pour ce cas)
field = sys.argv[2]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Début du calcul
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
id_fig = 1

for case in case_list:
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(case.case_name)
    # ----------------------------------------------------------------
    # Lecture du cas et information sur le nombre de cells rompues
    # ----------------------------------------------------------------
    hdf5_file = os.path.join(case.directory_name, "all_fields.hdf5")
    my_hd = OutputDatabaseExploit(hdf5_file)
    final_cell_status = my_hd.extract_field_at_time("CellStatus", my_hd.saved_times[-1])
    final_ruptured_cell_id = np.where(final_cell_status)[0]
    final_ruptured_cell_id = np.sort(final_ruptured_cell_id)

    # Vérification de la cohérence du champ demandé
    if field not in my_hd.saved_fields_type:
        print("The field type {:s} is not present in the database for case {:}. Available field types are {:s}!"
              .format(field, case.case_name, ",".join(my_hd.saved_fields_type)))

    # ----------------------------------------------------------------
    # Préparation de la figure
    # ----------------------------------------------------------------

    fig = plt.figure(id_fig)
    fig.suptitle(my_base_title[analysis_int], fontsize=20, fontweight='bold')
    # fig.suptitle(my_base_title[analysis_int] + " : " + case.case_name, fontsize=20, fontweight='bold')
    fig.patch.set_facecolor("white")
    n_colors = 500
    if analysis == "schlieren":
        plt.set_cmap('custom_cmap')  # pour utiliser la colormap qu'on a créé
    plt.xlabel(my_xlabel, fontsize=18)
    plt.ylabel(my_ylabel, fontsize=18)

    # ----------------------------------------------------------------
    # Construction de l'utility et construction de la map
    # ----------------------------------------------------------------
    utilities = SpaceTimeDiagramUtilities(case)
    X, Y, Z = utilities.build_XYZ_map_for_contourf_plot(field)


    vmin = np.min(Z)
    vmax = np.max(Z)
    if analysis == "schlieren":
        Z = np.gradient(Z, axis=0)
        # calcul de la dérivée temporelle (axis = 0) (un vecteur contenant tous les temps est Y[:,0]) de la pression Z
        vmin = -5.e8
        vmax = 5.e8
        print("Scaling the color map between {:} and {:}".format(vmin, vmax))

    # ----------------------------------------------------------------
    # Tracé color map 2D
    # ----------------------------------------------------------------
    if analysis in ["classic", "schlieren"] and plot:
        print("Plot the color map for analysis = {:}".format(analysis))

        if len(final_ruptured_cell_id) >= 1:
            # Séparation par paquets entre les mailles rompues pour prendre en compte le vide
            print("FYI : liste des cells rompues :" + str(final_ruptured_cell_id))
            ruptured_cell_id_after_offset = final_ruptured_cell_id + range(0, len(final_ruptured_cell_id))
            print("=> liste des cells rompues après décalage:" + str(ruptured_cell_id_after_offset))

            first_left_index = final_ruptured_cell_id[0]
            utilities.plot_section_of_color_map(X, Y, Z, end=first_left_index, n_colors=n_colors, vmin=vmin, vmax=vmax)
            offset = 1

            for i_rupture_index in ruptured_cell_id_after_offset[:-1]:
                right_current = i_rupture_index + 1
                left_next = ruptured_cell_id_after_offset[offset]
                utilities.plot_section_of_color_map(X, Y, Z, begin=right_current, end=left_next,
                                                    n_colors=n_colors, vmin=vmin, vmax=vmax, fig=id_fig)
                offset += 1

            last_right_index = ruptured_cell_id_after_offset[-1] + 1
            utilities.plot_section_of_color_map(X, Y, Z, begin=last_right_index, n_colors=n_colors,
                                                vmin=vmin, vmax=vmax, fig=id_fig)

        else:
            # Tracé tout simple
            utilities.plot_section_of_color_map(X, Y, Z, n_colors=n_colors, vmin=vmin, vmax=vmax, fig=id_fig)

        # Interface cible / projectile
        if utilities.data_has_interface():
            utilities.plot_interface(X, Y, 1)

        # Légende sous forme de colorbar
        ax, _ = matplotlib.colorbar.make_axes(fig.gca())
        colorbar = matplotlib.colorbar.ColorbarBase(ax, norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax))
        colorbar.set_label(my_colorbarlabel[analysis_int], fontsize=18)  # titre de la color bar
        # colorbar.set_label(r'\textit{bla}', fontsize=18)
        # colorbar.ax.set_title(field_name)

        id_fig += 1

        # fin de la boucle case pour analysis = classic ou schlieren

    # ----------------------------------------------------------------
    # Tracé d'une coupe :
    # ----------------------------------------------------------------
    if analysis == "coupe":

        # Redéfinition de quelques propriétés de la figure
        plt.figure(id_fig).suptitle("Coupe du diagramme space-time" + case.case_name, fontsize=20, fontweight='bold')
        plt.ylabel("{:}".format(field))

        # Coupe y = f(x) à t = constant
        if coupe_type == 0:
            plt.xlabel("Position [mm]")
            time = raw_input("Cas {:} :  t = ? (en secondes, "
                             "? pour afficher les temps existants dans les données) << " . format(case.case_name))
            # Pour afficher l'aide avant de reposer poliment la question
            if time == "?":
                print(my_hd.saved_times)
                time = raw_input("Cas {:} :  t = ? (en secondes, "
                                 "? pour afficher les temps existants dans les données) << ".format(case.case_name))
            time = float(time)
            index_time = np.where(np.array(my_hd.saved_times) <= time)[0][-1]
            # indice -1 pour prendre soit time soit le temps (inférieur) le plus proche
            plt.plot(X[index_time, :], Z[index_time, :], label=case.label, marker=case.marker, linestyle=case.linestyle)
            plt.title("Time : {:}".format(my_hd.saved_times[index_time]))

        # Coupe y = f(t) à item = constant (point de vue lagrangien)
        if coupe_type == 1:
            plt.xlabel("Temps [$\mu s$]")
            if len(final_ruptured_cell_id) > 0:
                ruptured_cell_id_after_offset = final_ruptured_cell_id + range(0, len(final_ruptured_cell_id))
                print("Cas {:} : les mailles rompues (après recalage) sont {:}".format(case.case_name,
                                                                                    ruptured_cell_id_after_offset))
            indice_item = int(raw_input("Cas {:} : item id = ? (entier après recalage au temps "
                                        "final) << ".format(case.case_name)))

            plt.plot(Y[:, indice_item], Z[:, indice_item], "-*", label=case.label + ' (cell ' + str(indice_item) + ')',
                     marker=case.marker, linestyle=case.linestyle)

    # ----------------------------------------------------------------
    # Affichage d'info sur un point à un temps donné :
    # ----------------------------------------------------------------
    if analysis == "pick":
        plot = False
        continue_bool = True

        if len(final_ruptured_cell_id) > 0:
            ruptured_cell_id_after_offset = final_ruptured_cell_id + range(0, len(final_ruptured_cell_id))
            print("Pour info, les mailles rompues (après recalage) sont " + str(ruptured_cell_id_after_offset))

        while continue_bool:
            time = raw_input("Cas {:} :  t = ? (en secondes) << ".format(case.case_name))
            time = float(time)
            index_time = np.where(np.array(my_hd.saved_times) <= time)[0][-1]

            indice_item = int(raw_input("Cas {:} : item id = ? (entier après recalage au temps "
                                        "final) << ".format(case.case_name)))

            print("----------------------------------")
            print("Item : {:}".format(indice_item))
            print("Temps : {:}".format(time))
            print("Position : {:}".format(Y[index_time, indice_item]))
            print("Valeur du champ {:} : {:}".format(field, Z[index_time, indice_item]))

            continue_bool = raw_input("Continuer ?") == "y"


if plot:
    plt.legend()
    plt.show()