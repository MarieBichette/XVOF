#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant une figure
"""
import matplotlib.pyplot as plt


class PhysicFigure(object):
    """
    Figure
    """
    def __init__(self, X, Y, xlabel="X", ylabel="Y", titre="titre"):
        self._fig = plt.figure()
        self._ax = self._fig.add_subplot(111)
        self._line, = self._ax.plot(X, Y, '-')
        self._ax.set_xlabel(xlabel)
        self._ax.set_ylabel(ylabel)
        self._ax.set_title(titre)
        plt.show(block=False)

    def set_y_limit(self, val_min=0., val_max=1.0):
        """ Fixation des limites en y"""
        self._ax.set_ylim([val_min, val_max])

    def update(self, X=None, Y=None):
        """ Mise à jour de l'image pour avoir une animation"""
        if(X is not None):
            self._line.set_xdata(X)
        if(Y is not None):
            self._line.set_ydata(Y)
        self._fig.canvas.draw()
        plt.show(block=False)
