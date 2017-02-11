#!/usr/bin/env python2.7
#  -*- coding: iso-8859-1 -*-
"""
Classe définissant les outils pour images temporelles

"""
import numpy as np
import os
import matplotlib.pyplot as plt


from xvof.figure_manager.time_figure import TimeFigure,Time_Field
from xvof.figure_manager.time_figure_tools import TimeFigureTools

class TimeFigureManager () :
    def __init__(self, simulation, item, id_item):
        '''
        :param simulation = nom de la simulation 'XFEM', 'REFERENCE' ou 'IMPOSEDPRESSURE'
        :param item : 'cell' or 'node'
        :param id_item: identifiant de cell or node of interest
        '''
        self.simulation = simulation
        self.item = item
        self.id_item = id_item


    def plot_single_time_figure(self):
        '''
        Plot a field history versus time
        '''
        time_fig_tools = TimeFigureTools(self.simulation, self.item, self.id_item)
        time_fig_tools.plot_data_history()


    def compare_time_fields(self):
        '''
        Plot the comparison between cell or node fields versus time
        '''
        ref_tools = TimeFigureTools("REFERENCE", self.item, self.id_item)
        pimp_tools = TimeFigureTools(self.simulation, self.item, self.id_item)

        ref_tools.plot_data_history()
        pimp_tools.plot_data_history()

    def compare_time_difference_fields(self):
        '''
        Plot the relative difference between ref and xfem cell or node fields versus time
        '''
        ref_tools = TimeFigureTools("REFERENCE", self.item, self.id_item)
        comp_tools = TimeFigureTools(self.simulation, self.item, self.id_item)

        ref_history = ref_tools.read_history_file()
        comp_history = comp_tools.read_history_file()

        ref_time = ref_history[:, 0]
        comp_time = comp_history[:, 0]
        verif_time = (ref_time == comp_time)

        if verif_time.all() == True:
            field_list = comp_tools.find_field_list()
            fig_num = 1
            for field in field_list:
                item_field = comp_tools.item_history_array[:, field.colonne_history]
                ref_item_field = ref_tools.item_history_array[:, field.colonne_history]
                # difference= np.abs((item_field - ref_item_field)/ref_item_field)*100.
                difference = np.abs(item_field - ref_item_field)
                ordre_de_grandeur_field = np.mean(np.abs(ref_item_field[np.where(ref_item_field > 1.)]))
                print 'Field {} ---------------------'.format(field.label)
                print 'Difference max = {}'.format(max(difference))
                print 'compared to reference {} ordre de grandeur field :{}'.format(field.label, ordre_de_grandeur_field)
                print '=> ratio {} % \n'.format(max(difference)/ordre_de_grandeur_field)
                item_figure = TimeFigure(ref_time, difference,self.id_item)
                item_figure.plot_time_figure(difference, field, comp_tools.set_graph_color(), comp_tools.simulation)

                plt.ylabel('Difference absolue avec REF')

        else:
            print "No comparison is possible because ref and simulation don't have the same length"






