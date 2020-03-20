#!/usr/bin/env python2.7
#  -*- coding: iso-8859-1 -*-
"""
Classe d�finissant les outils pour images temporelles

"""
import os

import matplotlib.pyplot as plt
import numpy as np

from xfv.src.figure_manager.time_figure import TimeFigure, Time_Field
from xfv.src.output_figure.hdf5_posttreatment_tools import read_database_in_array

# -- Rappel : Time_Field = namedtuple("Time_Field", ["label","title","application", "colonne_history"])
DensityField = Time_Field("Density [$kg/m^3$]", "Density", "cell", 1)
PressureField = Time_Field("Pressure [$Pa$]", "Pressure", "cell", 2)
EnergyField = Time_Field("Internal energy [$J/kg$]", "InternalEnergy", "cell", 3)
PseudoViscosityField = Time_Field("Pseudoviscosity [$Pa$]", "ArtificialViscosity", "cell", 4)
PositionField = Time_Field("Position [m]", "NodeCoordinates", "node", 5)
VelocityField = Time_Field("Vitesse [m/s]", "NodeVelocity", "node", 6)
DeviatoricStressField = Time_Field("Deviateur $S_{xx}$ [$Pa$]", "DeviatoricStress", "cell", 7)
StressField = Time_Field("contraintes $\sigma_{xx}$ [$Pa$]", "Stress", "cell", 8)
EquivalentPlasticStrainRate = Time_Field("Equivalent plastic strain rate", "EquivalentPlasticStrainRate", "cell", 9)
# les field.title correspondent aux noms des fields enregistr�s dans la bande hdf5

cell_field_list = [DensityField, PressureField, EnergyField, PseudoViscosityField,
                   DeviatoricStressField, StressField, EquivalentPlasticStrainRate]
node_field_list = [PositionField, VelocityField]
# ------------------------------------------


class TimeFigureTools:

    def __init__(self, case, item, id_number, data_file_name):
        """
        Define a class with methods to plot the history figures
        :param case : tuple Case
        :param item: 'cell' or 'node'
        :param id_number:  id of cell or node of interest
        """
        self.item = item
        self.id_number = id_number
        self.case = case
        self.history_directory = "./{:s}".format(case.directory_name)
        self._field_list = []
        self.data_file = data_file_name
        self.item_history_array = None

    def read_history_file(self, field):
        """
        Read the history file and extract an array with all history data
        field = PressureField par exemple
        :return: an array containing all time data
        """
        # import ipdb ; ipdb.set_trace()
        if self.data_file.split('.')[-1] == "dat":
            history_file = os.path.join(self.history_directory, self.data_file)
            print "Read data in {:}".format(history_file)
            reading = np.loadtxt(history_file, dtype='string', skiprows=2)
            time = reading[:, 0]
            field = reading[:, field.colonne_history]
            item_history = np.array(time, field)

        elif self.data_file.split('.')[-1] == "hdf5":
            field = field.title  # compatibilit� des deux noms de champs
            path_to_db = os.path.join(self.history_directory, self.data_file)
            item_history = read_database_in_array(path_to_db, self.id_number, field)
            print "Read data in {:}".format(path_to_db)

        else:
            raise ValueError(""" Wrong file extension specified """)
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

    def plot_data_history(self, field, my_label=None, my_color=None, subplot="1,1,1"):
        """
        plot all item fields versus time
        field : exemple PressureField
        :return:
        """
        if my_label is None:
            my_label = self.case.label
        if my_color is None:
            my_color = self.case.color
        self.read_history_file(field)
        time = self.item_history_array[:, 0]
        item_field_array = self.item_history_array[:, 1]
        item_figure = TimeFigure(time, item_field_array, self.id_number, subplot)
        item_figure.plot_time_figure(field, item_field_array, my_color, my_label, self.case.linestyle, self.case.marker)

    def save_history_figures(self, field, title_extension):
        """
        Save figures in the same repository as hdf5 band
        :param field : field of interest. ex: PressureField
        :param title_extension : str, precision on title of the figure + extension of save file(.png)
        """
        save_title = self.history_directory + field.title + "_{:}_{:}".format(self.item, self.id_number) + \
                     title_extension
        fig = plt.figure(field.colonne_history)
        fig.savefig(save_title)


def critere_couleur(critere, value):
    if value <= critere:
        color ='\033[92m'
    else:
        color = '\033[91m'
    return color
















