#!/usr/bin/env python2.7
# -*- coding: iso-8859-15 -*-
""" 
Compare the time evolution of nodes and cell for xfem simulation with reference
Parameters :
    USER INPUT :
    item : item to be studied (cell or node)
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

import matplotlib.pyplot as plt
import sys
from xvof.figure_manager.time_figure_manager import TimeFigureManager, set_acceptance_criterion, plot_error_vs_space
from xvof.figure_manager.time_figure_tools import TimeFigureTools
import xvof.case_definition as case_path

save_fig = True
show_fig = True
case_ref_choice = case_path.eps_025_ref_3x3
# case_ref_choice = case_path.default_ref
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

try:
    id_item = int(sys.argv[2])
except:
    raise ValueError("Le numéro d'item doit être un entier")

analysis = sys.argv[3]
if analysis not in ['single', 'gauche_droite', 'compare', 'diff', 'max_diff_vs_space']:
    raise ValueError("Type d'analysis {:s} inconnu! "
                     "Le type doit être parmi single, gauche_doite, compare, diff, max_diff_vs_space".format(analysis))

try:
    case_list = sys.argv[4].split(',')
except:
    case_list = [case_path.default]
    print "No case specified. Taking default database stored in path 0_XFEM / 0_REF\n" \
          "Enter a list of cases entre guillemet et sans espace"
# -----------------------------------------
# Initialize treatment and legend
# -----------------------------------------

if case_list == "all_mass":
    case_list = ["mass_1", "mass_2", "mass_3", "mass_4", "mass_5", "mass_6", "mass_7", "mass_8"]

if case_list == "special_mass":
    case_list = ["mass_1", "mass_3", "mass_5", "mass_7"]

if analysis == 'max_diff_vs_space':
    if id_item == -1:
        id_item_list = [500] # à modifier pour avoir quelque chose de régulier qui montre les erreurs
    else:
        id_item_list = [id_item]

plt.clf()
plt.close()

# -----------------------------------------
# Run user instruction (analysis)
# -----------------------------------------
for my_case in case_list:
    case = getattr(case_path, my_case)
    fig_manager = TimeFigureManager(case, item, id_item, case_ref=case_ref_choice)
    fig_tools = TimeFigureTools(case, item, id_item)

    if analysis == 'single':
        fig_manager.plot_single_time_figure()
        my_title_extension = "_{:}".format(getattr(case_path,case_list[0]).case_name) +".png"

    if analysis == 'gauche_droite':
        fig_manager.compare_gauche_droite()
        my_title_extension = "_gd.png"

    if analysis == 'compare':
        fig_manager.compare_time_fields()
        my_title_extension = ".png"

    if analysis == 'diff':
        print "Post treating data for case : {:}".format(case.label)
        fig_manager.plot_error_vs_time()
        my_title_extension = "_{:}_error.png".format(case_ref_choice.case_name)
        # set_acceptance_criterion(item, id_item, 0., 5., 1.e-4)

    if analysis == 'max_diff_vs_space':
        print "\033[31m//////// Case {:} //////////\033[37m".format(case.case_name)
        # Affichage dans le shell
        for id_item in id_item_list:
            xfem_manager = TimeFigureManager(case, item, id_item)
            print '\033[32m****** {} {} ******\033[37m'.format(xfem_manager.item, xfem_manager.id_item)
            xfem_manager.write_error_mass_case(critere=0.01)
        # Résultat graphique
        plot_error_vs_space(case, item, id_item_list, my_color=case.color, my_label=case.label)
        set_acceptance_criterion(item, id_item_list[0], min(id_item_list), max(id_item_list), 0.01)

if save_fig:
    fig_tools.save_history_figures(title_extension= my_title_extension)

if show_fig:
    plt.show()