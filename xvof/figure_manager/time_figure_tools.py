#!/usr/bin/env python2.7
#  -*- coding: iso-8859-1 -*-
"""
Classe définissant les outils pour images temporelles

"""
import numpy as np
import os


from xvof.figure_manager.time_figure import TimeFigure,Time_Field

#Rappel : Time_Field = namedtuple("Time_Field", ["label","title","application", "colonne_history"])
DensityField = Time_Field("Masse Volumique [kg/m3]", "Density","element", 1)
PressureField = Time_Field("Pression [Pa]", "Pressure","element", 2)
EnergyField = Time_Field("Energie interne [J/kg]", "Internal energy","element", 3)
PseudoViscosityField = Time_Field("Pseudoviscosite [Pa]", "Pseudo Viscosity","element", 4)

PositionField = Time_Field("Position [m]", "Position", "node",1)
VelocityField = Time_Field("Vitesse [m/s]", "Velocity","node", 2)

cell_field_list = [DensityField, PressureField, EnergyField, PseudoViscosityField]
node_field_list = [PositionField, VelocityField]

class TimeFigureTools() :

    def __init__(self, simulation, item, id_number):
        '''
        Define a class with methods to plot the history figures
        :param simulation: 'REFERENCE' or 'XFEM' or 'IMPOSED PRESSURE'
        :param item: 'cell' or 'node'
        :param id_number:  id of cell or node of interest
        '''
        self.simulation = simulation
        self.item = item
        self.id_number = id_number

    def read_history_file(self):
        '''
        Read the history file and extract an array with all history data
        :return: an array containing all time data
        '''
        history_file = ('./0_{:s}/{:s}_history_{:<03d}.dat').format(self.simulation.upper(), self.item, self.id_number)
        item_history = np.loadtxt(history_file, dtype='string', skiprows=3)
        item_history_array = item_history.astype(np.float)
        return item_history_array

    @property
    def item_history_array(self):
        return self.read_history_file()

    def set_graph_color(self):
        '''
        Choose a color to plot the graph versus the simulation considered
        :param simulation: REFERENCE , XFEM , IMPOSEDPRESSURE
        :return: color selected,str
        '''
        if self.simulation.upper() == "REFERENCE":
            sim_color = "red"
        elif self.simulation.upper() == "XFEM":
            sim_color = "blue"
        elif self.simulation.upper() == "IMPOSEDPRESSURE":
            sim_color = "green"
        else:
            print "Wrong simulation name : Possibilities are XFEM or IMPOSEDPRESSURE or REFERENCE !"
            sim_color = 'black'
        return sim_color

    def find_field_list(self):
        '''
        Define the field list for cell or node
        '''
        if self.item == "cell":
            field_list = cell_field_list
        elif self.item == "node":
            field_list = node_field_list
        else:
            print "Wrong type parameter in function plot_time_figure call." \
                  + os.linesep + "Only -cell- or -node- are possible"
        return field_list

    def plot_data_history(self):
        '''
        plot all item fields versus time
        :return:
        '''
        time = self.item_history_array[:,0]
        field_list = self.find_field_list()
        for field in field_list:
            item_field = self.item_history_array[:, field.colonne_history]
            item_figure = TimeFigure(time, item_field, self.id_number)
            item_figure.plot_time_figure(item_field, field, self.set_graph_color(), self.simulation )


