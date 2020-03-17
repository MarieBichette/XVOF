#!/usr/bin/env python2.7
#  -*- coding: iso-8859-1 -*-
"""
Classe définissant les outils pour images temporelles

"""
import matplotlib.pyplot as plt
import numpy as np
from pylab import rcParams

import xvof.utilities.case_definition as case_path
from xvof.figure_manager.time_figure import TimeFigure
from xvof.figure_manager.time_figure_tools import TimeFigureTools


class TimeFigureManager:
    def __init__(self, case, item, id_item, data_file_name, case_ref=None):
        """
        :param case : case to be considered
        :param item : 'cell' or 'node'
        :param id_item: identifiant de cell or node of interest
        :param data_file_name : name of the data results
        :param case_ref: case for reference plot
        """
        self.item = item
        self.id_item = id_item
        self.case = case
        self.case_ref = case_ref
        self.data_file = data_file_name
        # Création des case_tools
        self.time_fig_tools = TimeFigureTools(self.case, self.item, self.id_item, self.data_file)
        self.ref_tools = TimeFigureTools(self.case_ref, self.item, self.id_item, self.data_file)

    def plot_single_time_figure(self, field, my_label=None, my_color=None):
        """
        Plot a field history versus time
        """
        self.time_fig_tools.plot_data_history(field, my_label=my_label, my_color=my_color)

    # def compare_gauche_droite(self, field):
    #     """
    #     2 subplots in the figure to show the difference between cell 500 and 501
    #     """
    #     gaucheTools = TimeFigureTools(self.case, self.item, self.id_item -1, self.extension_file)
    #     gaucheTools.plot_data_history(field, my_label="Gauche", my_color="blue")
    #     # Partie droite déclarée dans le init de figure manager
    #     self.time_fig_tools.plot_data_history(field, my_label="Droite", my_color="orange")

    def compare_time_fields(self, field):
        """
        Plot the comparison between cell or node fields versus time
        """
        self.ref_tools.plot_data_history(field)
        self.time_fig_tools.plot_data_history(field)

    def compute_time_error_single(self, field, compute_error_index):
        """
        create an array with the max error (absolute difference between case and reference) in time for item
        cell or node number id_number_list and for case considered
        :return: array containing errors for cell or nodes fields for all times
        """
        self.ref_tools.read_history_file(field)
        self.time_fig_tools.read_history_file(field)
        item_field = self.time_fig_tools.item_history_array[:, 1]
        ref_item_field = self.ref_tools.item_history_array[:, 1]
        error_table = 0

        if compute_error_index == 0:
            print "Compute absolute error"
            error_table = np.abs(item_field - ref_item_field)

        if compute_error_index == 1:
            print "Compute erreur adimensionnée par max du champ (en valeur absolue)"
            error_table = (np.abs(item_field - ref_item_field)) / max(np.abs(ref_item_field - ref_item_field[0]))

        if compute_error_index == 2:
            print "Compute relative error"
            error_table = (np.abs(item_field - ref_item_field)) / np.abs(ref_item_field)

        return error_table

    # def write_error_mass_case(self, critere):
    #     """
    #     Ecrit dans le shell les résultats de l'erreur
    #     :param critere : float : critere de coloration des résultats (pour affichage en vert / rouge)
    #     """
    #     self.time_fig_tools.find_field_list()
    #     error = self.compute_time_error_single()
    #     for field in self.time_fig_tools.field_list:
    #         print '\033[34m------Field {} -----------------\033[37m'.format(field.label)
    #         max_error = max(error[:, field.colonne_history-1])
    #         couleur = critere_couleur(critere, max_error)
    #         print '=> erreur relative xfem / ref : {:}{:} \033[37m  \n'.format(couleur, max_error)

    def plot_error_vs_time(self, field, compute_error_index, my_label=None, my_color=None,
                           my_linestyle=None, my_marker=None):
        """
        Trace l'erreur entre xfem et ref en fonctiondu temps pour matrice masse = case
        :param field : field of interest
        :param compute_error_index : entier qui définit la méthode de calcul de l'erreur
        :param my_label : label for legend
        :param my_color : color for the graph
        :param my_linestyle: line style of the plot
        :param my_marker : markerof the data points
        """
        if my_color is None:
            my_color = self.case.color
        if my_label is None:
            my_label = self.case.label
        if my_linestyle is None:
            my_linestyle = self.case.linestyle
        if my_marker is None:
            my_marker = self.case.marker

        error = self.compute_time_error_single(field, compute_error_index)
        self.time_fig_tools.read_history_file(field)
        time = self.time_fig_tools.item_history_array[:, 0]

        item_figure = TimeFigure(time, error, self.id_item)
        item_figure.plot_time_figure(field, error, my_color, my_label, my_linestyle, my_marker)
        plt.ylabel('Erreur sur {:}'.format(field.title))


def enlarge_figure(inch_val_x, inch_val_y):
    rcParams['figure.figsize'] = inch_val_x, inch_val_y
