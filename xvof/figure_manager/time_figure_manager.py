#!/usr/bin/env python2.7
#  -*- coding: iso-8859-1 -*-
"""
Classe définissant les outils pour images temporelles

"""
import numpy as np
import matplotlib.pyplot as plt
from pylab import rcParams

from xvof.figure_manager.time_figure import TimeFigure
from xvof.figure_manager.time_figure_tools import TimeFigureTools, critere_couleur
import xvof.case_definition as case_path


class TimeFigureManager:
    def __init__(self, case, item, id_item, case_ref=case_path.default_ref):
        """
        :param case : case to be considered
        :param item : 'cell' or 'node'
        :param id_item: identifiant de cell or node of interest
        :param case_ref: case for reference plot
        """

        self.item = item
        self.id_item = id_item
        self.case = case
        self.simulation = self.case.simulation
        self.case_ref = case_ref

    def plot_single_time_figure(self):
        """
        Plot a field history versus time
        """
        time_fig_tools = TimeFigureTools(self.case, self.item, self.id_item)
        time_fig_tools.plot_data_history()

    def compare_gauche_droite(self):
        """
        2 subplots in the figure to show the difference between cell 500 and 501
        """
        gaucheTimeFigureTools = TimeFigureTools(self.case, self.item, self.id_item -1)
        gaucheTimeFigureTools.plot_data_history(my_label="Gauche", my_color="blue")
        droiteTimeFigureTools = TimeFigureTools(self.case, self.item, self.id_item)
        droiteTimeFigureTools.plot_data_history(my_label="Droite", my_color="orange")

    def compare_time_fields(self):
        """
        Plot the comparison between cell or node fields versus time
        """
        ref_tools = TimeFigureTools(self.case_ref, self.item, self.id_item)
        xfem_tools = TimeFigureTools(self.case, self.item, self.id_item)
        ref_tools.plot_data_history()
        xfem_tools.plot_data_history()

    def compute_time_error_adim_single(self):
        """
        create an array with the max error (absolute difference between case and reference) in time for item
        cell or node number id_number_list and for case considered
        :return: array containing errors for cell or nodes fields for all times
        """
        if self.id_item < 501:
            ref_tools = TimeFigureTools(self.case_ref, self.item, self.id_item)
        else:
            ref_tools = TimeFigureTools(self.case_ref, self.item, 500)
        xfem_tools = TimeFigureTools(self.case, self.item, self.id_item)
        xfem_tools.find_field_list()
        # Read file (va chercher les bons fichiers gràce à spécification de case dans TimeFigureTools
        ref_tools.read_history_file()
        xfem_tools.read_history_file()
        # Extract error for each field
        error_adim = np.zeros([np.shape(ref_tools.item_history_array)[0], np.shape(xfem_tools.field_list)[0]],
                              dtype=float)
        for field in xfem_tools.field_list:
            item_field = xfem_tools.item_history_array[:, field.colonne_history]
            ref_item_field = ref_tools.item_history_array[:, field.colonne_history]
            if self.id_item < 501:
                difference = np.abs(item_field - ref_item_field)
                ordre_de_grandeur_field = np.max(np.abs(ref_item_field))
                error_adim[:, field.colonne_history - 1] = difference / ordre_de_grandeur_field
            else:  # si moitié droite de la barre XFEM, on compare aux grandeurs initiales = init de cell 500 aussi
                difference = np.abs(item_field[200:])
                ordre_de_grandeur_field = np.max(np.abs(ref_item_field[200:]))
                error_adim[200:, field.colonne_history-1] = difference / ordre_de_grandeur_field
        return error_adim

    def write_error_mass_case(self, critere):
        """
        Ecrit dans le shell les résultats de l'erreur
        :param critere : float : critere de coloration des résultats (pour affichage en vert / rouge)
        """
        xfem_tools = TimeFigureTools(self.case, self.item, self.id_item)
        xfem_tools.find_field_list()
        error_adim = self.compute_time_error_adim_single()
        for field in xfem_tools.field_list:
            print '\033[34m------Field {} -----------------\033[37m'.format(field.label)
            max_error = max(error_adim[:, field.colonne_history-1])
            couleur = critere_couleur(critere, max_error)
            print '=> erreur adimensionnee par max de ref : {:}{:} \033[37m % \n'.format(couleur, max_error*100)

    def plot_error_vs_time(self, my_label=None, my_color=None):
        """
        Trace l'erreur entre xfem et ref en fonctiondu temps pour matrice masse = case
        :param my_label : label for legend
        :param my_color : color for the graph
        """
        if my_color is None:
            my_color = self.case.color
        if my_label is None:
            my_label = self.case.label
        xfem_tools = TimeFigureTools(self.case, self.item, self.id_item)
        xfem_tools.find_field_list()
        error_adim = self.compute_time_error_adim_single()
        xfem_tools.read_history_file()
        time = xfem_tools.item_history_array[:, 0]
        for field in xfem_tools.field_list:
            item_figure = TimeFigure(time, error_adim[:, field.colonne_history-1], self.id_item)
            item_figure.plot_time_figure(field, error_adim[:, field.colonne_history-1], my_color=my_color,
                                         my_label=my_label)
            plt.ylabel('Erreur adimensionnee par rapport a ref [-]')
            # plt.yscale('log')


def plot_error_vs_space(case, item, id_item_list, my_color, my_label):
    """
    Plot the error norm for each position and for defined mass case
    :param case: result case to be studied
    :param item: cell or node
    :param id_item_list: list of the item ids of interest
    :param my_color: color associated with the mass case
    :param my_label: label associated with the mass case
    """
    # initialisation de error_norm
    xfem_tools = TimeFigureTools(case, item, id_item_list[0])
    xfem_tools.find_field_list()
    error_norm = np.zeros([len(id_item_list), len(xfem_tools.field_list)])

    # remplissage de error_norm
    for id_item in id_item_list:
        xfem_manager = TimeFigureManager(case, item, id_item)
        error_adim_item = xfem_manager.compute_time_error_adim_single()
        for field in xfem_tools.field_list:
            error_norm[id_item_list.index(id_item), field.colonne_history-1] =\
                max(error_adim_item[:, field.colonne_history-1])

    # Tracé des résultats : error norm en fonction de item_id
    for field in xfem_tools.field_list:
        fig = plt.figure(field.colonne_history)
        plt.yscale('log')  # log scale
        fig.suptitle(field.title, fontsize=14, fontweight='bold')
        fig.patch.set_facecolor("white")
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(id_item_list, error_norm[:, field.colonne_history-1], color=my_color, label=my_label,
                linestyle='--', marker='*')
        ax.set_xlabel("Item id")  # à traduire en distance p/r discontinuité
        ax.set_ylabel("Erreur adimensionnee par rapport a ref [-]")


def set_acceptance_criterion(item, id_item, x_min, x_max, criterion_value):
    xfem_tools = TimeFigureTools("mass_1", item, id_item)  # juste pour avoir la taille des champs à tracer
    xfem_tools.find_field_list()
    for id_fig in range(1, len(xfem_tools.field_list)+1):
        fig = plt.figure(id_fig)
        ax = fig.add_subplot(1, 1, 1)
        ax.plot([x_min, x_max], [criterion_value, criterion_value],
                color="red", label="CRITERE", linestyle='-', linewidth=3)
        ax.legend(loc='best')
        # Pour avoir une légende sur le coté
        # box = ax.get_position()
        # ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
        # ax.legend(loc='center left', bbox_to_anchor=(1., 0.5))


def enlarge_figure(inch_val_x, inch_val_y):
    rcParams['figure.figsize'] = inch_val_x, inch_val_y
