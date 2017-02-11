#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant une figure en fonction du temps
"""
import matplotlib.pyplot as plt


from collections import namedtuple

from xvof.data.data_container import DataContainer


Time_Field = namedtuple("Time_Field", ["label","title","application", "colonne_history"])

class TimeFigure():
    """
    Figure (en fonction du temps)
    = cree la figure à partir de la classe mère mais avec le temps en abscisse
    Y : champs à tracer
    """
    def __init__(self, temps=[], champ_Y=[], id_number=-1):
        self._temps = temps
        self._champ_Y = champ_Y
        self._id_number = id_number


    def plot_time_figure(self, field_Y, field_name, my_color, my_sim_label):
        """
        plot the figure vecteur_Y = f(t)
        :var field_Y : champ à tracer
                    = array contenant les valeurs d'un champ
        :var field_name : nom du champ à tracer (de type Time_Field)
        :var item_id : numéro de l'item (cell ou node) à tracer
        :var mycolor, sim_label : paramètres couleur et légende pour la figure
        """
        self._temps = 1.e+06 * self._temps

        plt.figure(field_name.colonne_history)
        plt.plot(self._temps, field_Y, color = my_color, label = my_sim_label, marker = '.',linestyle='-' )
        plt.xlabel("time [micro s]")
        plt.ylabel(field_name.label)
        plt.title(field_name.title + " in " + field_name.application + " "+ str(self._id_number))
        plt.legend(bbox_to_anchor=(1,1), bbox_transform=plt.gcf().transFigure)

