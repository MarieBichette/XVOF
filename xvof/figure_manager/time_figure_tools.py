#!/usr/bin/env python2.7
#  -*- coding: iso-8859-1 -*-
"""
Classe définissant les outils pour images temporelles

"""
import numpy as np
import os
import matplotlib.pyplot as plt

from xvof.figure_manager.time_figure import TimeFigure, Time_Field

# -- Rappel : Time_Field = namedtuple("Time_Field", ["label","title","application", "colonne_history","min","max"])
DensityField = Time_Field("Density [$kg/m^3$]", "Density","cell", 1, 8124, 8130)
PressureField = Time_Field("Pressure [$Pa$]", "Pressure","cell", 2, -7.e+07, 2.e+07)
EnergyField = Time_Field("Internal energy [$J/kg$]", "InternalEnergy","cell", 3, 7.5, 10.)
PseudoViscosityField = Time_Field("Pseudoviscosity [$Pa$]", "ArtificialViscosity","cell", 4, 0., 700000.)
# PositionField = Time_Field("Position [m]", "NodeCoordinates", "node",1, 0.0120, 0.0125)
VelocityField = Time_Field("Vitesse [m/s]", "NodeVelocity","node", 1, -700., 100.)

cell_field_list = [DensityField, PressureField, EnergyField, PseudoViscosityField]
node_field_list = [VelocityField]
# ------------------------------------------


class TimeFigureTools:

    def __init__(self, case, item, id_number):
        """
        Define a class with methods to plot the history figures
        :param case : tuple Case
        :param item: 'cell' or 'node'
        :param id_number:  id of cell or node of interest
        """
        self.item = item
        self.id_number = id_number
        self.case = case
        self.simulation = self.case.simulation
        self.history_directory = "./0_{:s}/{:s}".format(self.simulation.upper(), case.directory_name)
        self._field_list = []

    def read_history_file(self):
        """
        Read the history file and extract an array with all history data
        :return: an array containing all time data
        """
        data_file = "{:s}_history_{:<03d}.dat".format(self.item, self.id_number)
        history_file = os.path.join(self.history_directory, data_file)
        item_history = np.loadtxt(history_file, dtype='string', skiprows=2)
        self.item_history_array = item_history.astype(np.float)

    def find_field_list(self):
        """
        Define the field list for cell or node
        """
        if self.item == "cell":
            self._field_list = cell_field_list
        elif self.item == "node":
            self._field_list = node_field_list
        else:
            print "Wrong type parameter in function plot_time_figure call." \
                  + os.linesep + "Only -cell- or -node- are possible"
        return self._field_list

    @property
    def field_list(self):
        return self._field_list

    def plot_data_history(self, my_label=None, my_color=None, subplot="1,1,1"):
        """
        plot all item fields versus time
        :return:
        """
        if my_label is None:
            my_label = self.case.label
        if my_color is None:
            my_color = self.case.color
        self.read_history_file()
        time = self.item_history_array[:, 0]
        self.find_field_list()
        for field in self.field_list:
            item_field_array = self.item_history_array[:, field.colonne_history]
            item_figure = TimeFigure(time, item_field_array, self.id_number, subplot)
            item_figure.plot_time_figure(field, item_field_array, my_color, my_label)

    def save_history_figures(self, title_extension):
        """
        Save figures in the same repository as hdf5 band
        :param title_extension : str, precision on title of the figure + extension of save file(.png)
        """
        self.find_field_list()
        for field in self.field_list:
            save_title = self.history_directory + field.title + "_{:}_{:}".format(self.item, self.id_number) + \
                         title_extension
            fig = plt.figure(field.colonne_history)
            fig.savefig(save_title)


def critere_couleur(critere, value):
    if value <= critere:
        couleur='\033[92m'
    else:
        couleur ='\033[91m'
    return couleur
















