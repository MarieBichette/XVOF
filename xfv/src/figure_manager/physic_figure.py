#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe d�finissant une figure
"""
import matplotlib.pyplot as plt
from os import sep
import numpy as np
from xfv.src.data.data_container import DataContainer


class PhysicFigure(object):
    """
    Figure
    """
    def __init__(self, X, Y, xlabel="X", ylabel="Y", titre="titre", interface_id=0, save_path=None):
        self._fig, self._ax = plt.subplots()
        self._line_target, = self._ax.plot(X, Y, '-+', color='blue')
        self._line_projectile, = self._ax.plot(X, Y, '-+', color='purple')
        self._ax.set_xlabel(xlabel)
        self._ax.set_ylabel(ylabel)
        self._ax.set_title(titre)
        self._fig_number = 1
        self._title = titre
        self._fig.canvas.draw()
        self._save_path = save_path
        self.interface_id = interface_id # pour s�parer cible et projectile
        plt.show(block=False)

        # x = np.array([0, 0.01])
        # y = np.array([DataContainer().material_target.failure_model.failure_criterion_value,
        #               DataContainer().material_target.failure_model.failure_criterion_value])
        # self._line_ref, = self._ax.plot(x, y, color='red')


    def set_y_limit(self, val_min=0., val_max=1.0):
        """ Fixation des limites en y"""
        self._ax.set_ylim([val_min, val_max])

    def set_x_limit(self, val_min=0., val_max=1.0):
        """ Fixation des limites en x"""
        self._ax.set_xlim([val_min, val_max])

    def update(self, X=None, Y=None, title_comp=None):
        """
        Mise � jour de l'image pour avoir une animation
        Sauvegarde de l'image si pr�sence d'un path
        """
        if X is not None:
            self._line_target.set_xdata(X[self.interface_id:])
            self._line_projectile.set_xdata(X[:self.interface_id])
        if Y is not None:
            self._line_target.set_ydata(Y[self.interface_id:])
            self._line_projectile.set_ydata(Y[:self.interface_id])
        if title_comp is not None:
            self._ax.set_title(self._title + ' ' + title_comp)

        self._ax.draw_artist(self._ax.patch)
        self._ax.draw_artist(self._line_target)
        self._ax.draw_artist(self._line_projectile)
        # self._ax.draw_artist(self._line_ref)
        self._fig.canvas.update()
        self._fig.canvas.flush_events()

        if self._save_path is not None:
            rac_path = self._save_path + sep + self._title
            fig_path = rac_path + "_{:04d}.png".format(self._fig_number)
            fig_path = fig_path.replace(" ", "_")
            data_path = rac_path + "_{:04d}.dat".format(self._fig_number)
            data_path = data_path.replace(" ", "_")
            self._fig.savefig(fig_path)
            self._fig_number += 1
            with open(data_path, "w") as fo:
                for a,b in zip(X, Y):
                    fo.write("{:20.18g}{:s}{:20.18g}\n".format(float(a), " ", float(b)))


