#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant une figure en fonction du temps
"""
import matplotlib.pyplot as plt
from collections import namedtuple


Time_Field = namedtuple("Time_Field", ["label", "title", "application", "colonne_history"])


class TimeFigure:
    """
    Figure (en fonction du temps)
    = cree la figure à partir de la classe mère mais avec le temps en abscisse
    Y : champs à tracer
    """
    def __init__(self, temps, champ_y, id_number=-1, subplot="1,1,1"):
        self._temps = temps
        self._champ_Y = champ_y
        self._id_number = id_number
        self._subplot = subplot.split(',')

    def plot_time_figure(self, field_name, field_y, my_color, my_label, my_linestyle, my_marker):
        """
        plot the figure vecteur_Y = f(t)
        :var field_y : champ à tracer
                    = array contenant les valeurs d'un champ
        :var field_name : nom du champ à tracer (de type Time_Field)
        :var my_color, my_label : paramètres couleur et légende pour la figure
        """
        fig = plt.figure(field_name.colonne_history, figsize=(9,7))
        fig.suptitle(field_name.title, fontsize=20, fontweight='bold')
        fig.patch.set_facecolor("white")
        ax = fig.add_subplot(1,1,1)
        ax.plot(self._temps * 1.e+06, field_y, color=my_color, label=my_label, linestyle=my_linestyle, marker=my_marker)
        ax.set_xlabel("Time [$\mu s$]", fontsize=20)
        ax.set_title(field_name.application + " " + str(self._id_number))

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 1., box.height])

        ax.set_ylabel(field_name.label, fontsize=20)
        # ax.set_ylim([8920, 8932])
        # ax.set_xlim([2.5, 5.])
        # ax.legend(loc=4, fontsize=15)
        # plt.legend(bbox_to_anchor=(0.5, 1), loc=2,  mode="expand", borderaxespad=0.)
