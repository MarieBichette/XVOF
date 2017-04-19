#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant une figure
"""
import matplotlib
from os import sep


class PhysicFigure(object):
    """
    Figure
    """
    def __init__(self, X, Y, xlabel="X", ylabel="Y", titre="titre"):
        self._fig, self._ax = matplotlib.pyplot.subplots()
        self._line, = self._ax.plot(X, Y, '-+')
        self._ax.set_xlabel(xlabel)
        self._ax.set_ylabel(ylabel)
        self._ax.set_title(titre)
        self._fig_number = 1
        self._title = titre
        self._fig.canvas.draw()
        matplotlib.pyplot.show(block=False)

    def set_y_limit(self, val_min=0., val_max=1.0):
        """ Fixation des limites en y"""
        self._ax.set_ylim([val_min, val_max])

    def set_x_limit(self, val_min=0., val_max=1.0):
        """ Fixation des limites en x"""
        self._ax.set_xlim([val_min, val_max])

    def update(self, X=None, Y=None, title_comp=None):
        """
        Mise à jour de l'image pour avoir une animation
        Sauvegarde de l'image si présence d'un path
        """
        if(X is not None):
            self._line.set_xdata(X)
        if(Y is not None):
            self._line.set_ydata(Y)
        if(title_comp is not None):
            self._ax.set_title(self._title + ' ' + title_comp)

        self._ax.draw_artist(self._ax.patch)
        self._ax.draw_artist(self._line)
        self._fig.canvas.update()
        self._fig.canvas.flush_events()

