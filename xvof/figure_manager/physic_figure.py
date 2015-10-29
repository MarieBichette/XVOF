#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant une figure
"""
from os import sep

import matplotlib.pyplot as plt


class PhysicFigure(object):
    """
    Figure
    """
    def __init__(self, X, Y, xlabel="X", ylabel="Y", titre="titre", save_path=None):
        self._fig = plt.figure()
        self._ax = self._fig.add_subplot(111)
        self._line, = self._ax.plot(X, Y, '-+')
        self._ax.set_xlabel(xlabel)
        self._ax.set_ylabel(ylabel)
        self._ax.set_title(titre)
        self._save_path = save_path
        self._fig_number = 1
        self._title = titre

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
        self._fig.canvas.draw()
        if (self._save_path is not None):
            rac_path = self._save_path + sep + self._title
            fig_path = rac_path + "_{:04d}.png".format(self._fig_number)
            fig_path = fig_path.replace(" ", "_")
            data_path = rac_path + "_{:04d}.dat".format(self._fig_number)
            data_path = data_path.replace(" ", "_")
            self._fig.savefig(fig_path)
            self._fig_number += 1
            with open(data_path, 'w') as fo:
                for a, b in zip(X, Y):
                    fo.write("{:20.18g}{:s}{:20.18g}\n".format(float(a), 4 * " ", float(b)))

