#!/usr/bin/env python2.7
# -*- coding: utf-8-*-

"""
Script pour tracer les vitesses de surface libre étude paramtétrique du
modèle cohésif avec les données expérimentales TA5
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from xvof.output_figure.hdf5_posttreatment_tools import read_database_in_array
from xvof.utilities.experimental_data_exploit import ExperimentalData

msg = "Petit programme tout mignon pour étude param deslois cohésives \n"
msg += "Ce script attend comme  arguments  : \n"
msg += "- un paramètre à regarder : strength, separation, slope\n"
msg += "- -h ou --help pour afficher l'aide\n"

if sys.argv[1] in ["-h", "--help"]:
    print msg
    exit(0)

if len(sys.argv) != 2:
    print msg
    exit(0)

analysis = sys.argv[1]

def plot_cohesive_law(id_fig, figure_subtitle, cohesive_parameters, color_list):
    """
    Trace la loi cohésive utilisée
    :param id_fig: int, figure id which to plot in
    :param figure_subtitle: str, subtitle for the figure
    :param cohesive_parameters: array [n*3] avec les 3 paramètres des lois cohésivesà tracer
    :return:
    """
    fig = plt.figure(id_fig)
    fig.patch.set_facecolor("white")
    fig.suptitle("Cohesive laws $\sigma = f(\delta)$", fontsize=20, fontweight='bold')
    plt.title(figure_subtitle)
    plt.xlabel("Discontinuity opening $\delta $ [mm]", fontsize=18)
    plt.ylabel("Cohesive stress $\sigma$ [GPa]", fontsize=18)

    for i in range(0, np.shape(cohesive_parameters)[0]):
        sigma = cohesive_parameters[i, 0]
        crit = cohesive_parameters[i, 1] * 1000  # pour passer en mm
        plat = cohesive_parameters[i, 2] * 1000  # pour passer en mm
        label = "$\sigma$ {:} crit {:} plateau {:}".format(sigma, crit, plat)
        plt.plot([0, plat, crit], [sigma, sigma, 0], linewidth=2, label=label, color=color_list[i])
        plt.ylim(0, 1.1 * sigma)


def plot_series_of_data(id_fig, figure_subtitle, exp_data_path, array_of_data, color_list):
    """
    Read and the hd5f data and plot the free surface velocity
    :param id_fig : int, figure id which to plot in
    :param figure_subtitle: str, subtitle for the figure
    :param exp_data_path: path to the experimental data
    :param array_of_data: list of data
    [[sigma, delta_crit, delta_plateau],[...]] to simulation hdf5 results
    :return:
    """
    # Préparation de la figure
    fig = plt.figure(id_fig)
    fig.patch.set_facecolor("white")
    # fig.suptitle("Vitesse de surface libre TA 5", fontsize=20, fontweight='bold')
    fig.suptitle("Free surface velocity", fontsize=20, fontweight='bold')
    plt.title(figure_subtitle)
    plt.xlabel("Time [$\mu$s]", fontsize=18)
    plt.ylabel("Velocity of free surface [m/s]",fontsize=18)

    # Tracé des courbes expérimentales
    exp_result = ExperimentalData(exp_data_path)
    exp_result.plot_experimental_results(id_fig)

    # Tracé des résultats de la simulation
    for i in range(0, np.shape(array_of_data)[0]):
        sigma = array_of_data[i, 0]
        crit = array_of_data[i, 1]
        plat = array_of_data[i, 2]
        path = \
            "./ETUDE_PARAMETRIQUE_TA5/Cohesive_strength_{:}MPa/Critical_separation_{:.0E}m/Plateau_at_{:.1E}".\
            format(int(sigma * 1000), crit, plat)
        file = os.path.join(path, "all_fields.hdf5")
        print "Reading data in {:}".format(file)
        hdf5_data = read_database_in_array(file, -1, "NodeVelocity")
        time = hdf5_data[:, 0]
        velocity = hdf5_data[:, 1]
        mask = (velocity > 1)
        time_0 = time[np.where(mask)[0][0]]  # recalage en temps pour le TA5
        label = "$\sigma$ {:} crit {:} plateau {:}".format(sigma, crit, plat)
        plt.plot((time - time_0) * 1.e+06, velocity, label=label, color=color_list[i])
        # plt.plot((time - time_0) * 1.e+06, velocity, label=label, color=color_list[i])
    # plt.legend(loc="best")


experimental_file = "TA5.TXT"
experimental_dir = "//home/marie/Documents/These/Experimental_data/"
exp_data = os.path.join(experimental_dir, experimental_file)

# for sigma in sigma_s:
#     print "-------- Analysis for sigma max = {:} GPa".format(sigma)
#     for crit in critic:
#         print "--- Analysis for critical opening = {:} ".format(crit)
#         for plat in plateau:
#             print "- Analysis for plateau = {:} ".format(plat)

color_list= ["blue", "red", "orange", "green"]
# Influence de la cohesive strength
if analysis == "strength":
    print "Post-treatment : cohesive strength"
    array_strength_1 = np.array([[1, 1.e-4, 0.8e-4], [2, 1.e-4, 0.6e-4], [3, 1.e-4, 0.4e-4]])
    array_strength_2 = np.array([[1, 1.e-4, 0.9e-4], [2, 1.e-4, 0.8e-4], [3, 1.e-4, 0.7e-4]])
    array_strength_3 = np.array([[1, 1.e-4, 1.e-4], [2, 1.e-4, 1.e-4], [3, 1.e-4, 1.e-4]])
    # plot_series_of_data(1, "Influence of the cohesive strength (pente faible)", exp_data, array_strength_1, color_list)
    # plot_series_of_data(2, "Influence of the cohesive strength (pente elevee)", exp_data, array_strength_2, color_list)
    plot_series_of_data(2, " ", exp_data, array_strength_2, color_list)
    # plot_series_of_data(3, "Influence of the cohesive strength (pente infinie)", exp_data, array_strength_3, color_list)

    # plot_cohesive_law(11, "Influence of the cohesive strength (pente faible)", array_strength_1, color_list)
    plot_cohesive_law(12, "Influence of the cohesive strength (pente elevee)", array_strength_2, color_list)
    plot_cohesive_law(12, " ", array_strength_2, color_list)
    # plot_cohesive_law(13, "Influence of the cohesive strength (pente infinie)", array_strength_3, color_list)

# Influence de la critical separation
if analysis == "separation":
    print "Post-treatment : critical separation"
    array_separation_1 = np.array([[2, 0.5e-4, 0.1e-4], [2, 0.6e-4, 0.2e-4], [2, 0.8e-4, 0.4e-4], [2, 1.e-4, 0.6e-4]])
    array_separation_2 = np.array([[2, 0.5e-4, 0.3e-4], [2, 0.6e-4, 0.4e-4], [2, 0.8e-4, 0.6e-4], [2, 1.e-4, 0.8e-4]])
    array_separation_3 = np.array([[2, 0.5e-4, 0.5e-4], [2, 0.6e-4, 0.6e-4], [2, 0.8e-4, 0.8e-4], [2, 1.e-4, 1.e-4]])
    # plot_series_of_data(4, "Influence of the critical separation (pente faible)", exp_data, array_separation_1, color_list)
    # plot_series_of_data(5, "Influence of the critical separation (pente elevee)", exp_data, array_separation_2, color_list)
    plot_series_of_data(5, " ", exp_data, array_separation_2, color_list)
    # plot_series_of_data(6, "Influence of the critical separation (pente infinie)", exp_data, array_separation_3, color_list)

    # plot_cohesive_law(14, "Influence of the critical separation (pente faible)", array_separation_1, color_list)
    # plot_cohesive_law(15, "Influence of the critical separation (pente elevee)", array_separation_2, color_list)
    plot_cohesive_law(15, " ", array_separation_2, color_list)
    # plot_cohesive_law(16, "Influence of the critical separation (pente infinie)", array_separation_3, color_list)

if analysis == "slope":
    print "Post-treatment : slope (with constant dissipated energy)"
    array_slope_1 = np.array([[1, 0.6e-4, 0.6e-4], [1, 0.8e-4, 0.4e-4], [1, 1.e-4, 0.2e-4]])
    array_slope_2 = np.array([[2, 0.6e-4, 0.6e-4], [2, 0.8e-4, 0.4e-4], [2, 1.e-4, 0.2e-4]])
    array_slope_3 = np.array([[3, 0.6e-4, 0.6e-4], [3, 0.8e-4, 0.4e-4], [3, 1.e-4, 0.2e-4]])
    plot_series_of_data(7, "Influence of the softening slope (stress = 1)", exp_data, array_slope_1, color_list)
    plot_series_of_data(8, "Influence of the softening slope (stress = 2)", exp_data, array_slope_2, color_list)
    plot_series_of_data(9, "Influence of the softening slope (stress = 3)", exp_data, array_slope_3, color_list)

    plot_cohesive_law(17, "Influence of the softening slope (stress = 1)", array_slope_1, color_list)
    plot_cohesive_law(18, "Influence of the softening slope (stress = 2)", array_slope_2, color_list)
    plot_cohesive_law(19, "Influence of the softening slope (stress = 3)", array_slope_3, color_list)

#
plt.show()



