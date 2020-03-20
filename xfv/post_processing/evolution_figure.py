#!/usr/bin/env python2.7
# -*- coding: iso-8859-15 -*-
""" 
Compare the time evolution of nodes and cell for xfem simulation with reference
Parameters :
    USER INPUT :
    id_item : id of the item to be studied (int)
        For 'max_diff_vs_space' : possibility to consider a list of ids : id_item = -1
            -> takes into account a list entered in the source code
    analysis : type of analysis to be performed :
            - single : plot the item fields for item item n° id_item
            - gauche_droite : superpose the left and right part of the ruptured cell
                    (other parameters are fixed : id_item = 501)
            -compare : superpose the case simulation item fields and the reference fields
                    (default reference case = masse_wilkins)
            - diff : plot the error vs time for case case (error calculated with ref_case)
            - max_diff_vs_space : plot max error in time (computed with ref_case) for all id_item in id_item_list
    case_list : list of cases to be considered.
            case : tuple qui contient : nom du cas, simulation associée,
                    répertoire de stockage des résultats, label et couleur pour les graphiques
            entrer cette liste entre guillemet et sans espace. séparateur = ','.
            case default : pour prendre les résultats situés dasn "0_XFEM"

    OTHER PARAMETERS
    save_fig, show_fig : booléens pour sauver et afficher les résultats
    id_item_list : cf id_item
    ref_case : default : mass_wilkins
"""

import os
import sys

import matplotlib.pyplot as plt
from matplotlib2tikz import save as tikz_save

import xfv.src.figure_manager.time_figure_tools as fig_tools_path
from xfv.src.utilities.case_definition import CaseManager
from xfv.src.figure_manager.time_figure_manager import TimeFigureManager
from xfv.src.figure_manager.time_figure_tools import TimeFigureTools
from xfv.src.output_figure.profile_tools import A3_list, plot_field_from_txt_file, get_field_from_txt_file

# -----------------------------------------
# Initialize treatment and legend
# -----------------------------------------
save_fig = True
show_fig = True
file_name = "all_fields.hdf5"
extension = '.hdf5'

project_dir = os.path.split(os.path.dirname(os.path.abspath(os.path.curdir)))[0]
base_dir = os.path.split(project_dir)[0]

# -----------------------------------------
# Read user instructions
# -----------------------------------------
msg = "Mini programme pour tracer l'évolution temporelle d'un champ en un point de la barre\n" \
      "Ce script attend comme arguments  : \n"
msg += "- le type de simulation (single|compare|diff|compare_with_A3) \n"
msg += "- les champs à tracer ([X]Field, ...)\n"
msg += "- les ids des items (int, séparés par une virgule et sans espace) \n"
msg += "- une liste des cas à traiter (séparés par une virgule, sans espace | None par défaut) \n"
msg += "- -h ou --help pour afficher l'aide\n"
msg += " Si analysis = diff : " \
    "\n \t - Le premier cas rensigné sert de référence" \
    "\n \t - On demandera une option pour savoir comment calculer l'erreur (absolue, adimensionnée, relative)"

analysis = None

if len(sys.argv) not in [2, 5]:
    print msg
    exit(0)

if sys.argv[1] in ["-h", "--help"]:
    print msg
    exit(0)

else: # vraie simulation : on lit lest arguments données par l'utilisateur
    try:
        analysis = sys.argv[1]
        if analysis not in ['single', 'compare', 'diff', 'compare_with_A3']:
            print msg
            print("Type d'analysis {:s} inconnu!".format(analysis))
            os._exit(0)

        field_list = sys.argv[2].split(',')
        field_list = [getattr(fig_tools_path, t_f) for t_f in field_list]

        id_item_list = sys.argv[3].split(',')
        id_item_list = [int(id_item) for id_item in id_item_list]

    except:
        print msg
        raise

case_list = sys.argv[4].split(',')
case_list = [CaseManager().find_case(case) for case in case_list]
case_ref_choice = case_list[0]  # pas utile pour une comparaison single mais on a beosin d'en définir un quand même

if analysis in ['diff']:
    # Pour une courbe d'erreur, la référence est donnée par le premier cas renseigné
    if len(case_list) < 2:
        print "Impossible de faire une analyse en erreur (diff) si il n'y a pas au moins 2 cas renseignés"
        exit(0)
    case_ref_choice = case_list[0]
    case_list = case_list[1:]
    # On demande quel calcul d'erreur faire
    option = 1
    # option = None
    # while option not in range(3):
    #     option = input("Quelle méthode de calcul pour l'erreur ? (0 = absolue, 1 = adimensionnée, 2 = relative) ")

if analysis=="compare_with_A3":
    # uniquement certaines données sont actuellement disponibles pour la comparaison avec A3
    if id_item not in [0, 100, 200, 500, 1000]:
        raise ValueError("Data not available for cell {:}\n \
        Enter a cell in [0, 100, 200, 500, 1000]".format(id_item))

# -----------------------------------------
# Run user instruction (analysis)
# -----------------------------------------
for case in case_list:
    print "-------------------------------"
    print case.case_name
    item = field_list[0].application

    for field in field_list:
        print field.title

        for id_item in id_item_list:
            print "item : {:}".format(id_item)

            fig_manager = TimeFigureManager(case, item, id_item, data_file_name=file_name, case_ref=case_ref_choice)
            fig_tools = TimeFigureTools(case, item, id_item, data_file_name=file_name)

            if analysis == 'single':
                fig_manager.plot_single_time_figure(field)
                plt.legend(loc="best")
                if save_fig:
                    tikz_save(os.path.join(base_dir,
                                           "Documents/These/Matplotlib2Tikz/time_evolution_{:}_{:}_{:}.tex".format(field.title,
                                                                                                            item,
                                                                                                            id_item)))

            elif analysis == 'compare':
                fig_manager.compare_time_fields(field)
                plt.legend(loc=4)

            elif analysis == 'diff':
                print "Reference case is {}".format(case_ref_choice.case_name)
                print "Post treating data for case : {:}".format(case.label)
                fig_manager.plot_error_vs_time(field, compute_error_index=option)
                # set_acceptance_criterion(item, id_item, 0., 5., 1.e-4)
                plt.legend(loc=4)

                if save_fig:
                    tikz_save(os.path.join(base_dir,
                                           "Documents/These/Matplotlib2Tikz/time_evolution_{:}_{:}_{:}.tex".format(field.title,
                                                                                                            item,
                                                                                                            id_item)))

            elif analysis == 'compare_with_A3':
                # ^^^^^^ Traitement de la bande hdf5 :^^^^^^^^
                fig_manager.plot_single_time_figure(field)
                # ^^^^^^ Ajout des mesures A3 ^^^^^^^^
                directory = "//home/marie/Documents/These/export_A3/"  # Find directory of A3 data
                subdir = {'elasto': 'RESULT_ELASTO/', 'epp_120_MPa': 'RESULT_EPP_120_MPa/',
                          'hydro_sans_pseudo': 'RESULT_HYDRO_sans_pseudo', 'hydro_avec_pseudo': 'RESULT_HYDRO'}
                directory = os.path.join(directory, subdir[case], 'evolution')
                path_to_file = os.path.join(directory,
                                            "evolution_{:s}_{}.txt".format(A3_list[field.title], id_item))
                x_a3, y_a3 = get_field_from_txt_file(path_to_file)
                plot_field_from_txt_file(field.colonne_history, x_a3, y_a3, multiplicateur=1.e+06)
                plt.xlabel("Time [$\mu s$]")
                plt.ylabel(field.label)
                plt.legend(loc=4)

if show_fig:
    plt.show()




