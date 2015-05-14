#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe dÃ©finissant le gestionnaire d'images

@todo: Hériter de la classe Figure
@todo: Transformer la classe Field en NamedTuple (cf. properties)
"""
from collections import namedtuple
from os import makedirs
from os.path import exists

import matplotlib.pyplot as plt
import numpy as np
from physic_figure import PhysicFigure


Field = namedtuple("Field", ["label", "titre", "val_min", "val_max", "results_path"])

PressureField = Field("Pression [Pa]", "Champ de pression", -7.5e+09, 7.5e+09, "./RESULTATS/PressureField")
DensityField = Field("Masse volumique [kg/m3]", "Champ de densite", 7500.0, 8500.0, "./RESULTATS/DensityField")
InternalEnergyField = Field("Energie interne [J/kg]", "Champ d energie interne", 0, 40000.0, "./RESULTATS/InternalEnergyField")
PseudoViscosityField = Field("Pseudoviscosite [Pa]", "Champ de pseudoviscosite", 0, 1.0e+09, "./RESULTATS/PseudoViscosityField")
CellPositionField = Field("Position [m]", "Champ de position", 0, 0.02, "./RESULTATS/CellPositionField")

class FigureManager(object):
    """
    Gestionnaire de figures
    """
    def __init__(self, mesh_instance, dump=False, show=True):
        self.__mesh_instance = mesh_instance
        self.__figures_mailles = []
        self.__champs_mailles = {}
        self.update_fields()
        self._dump = dump
        self._show = show

    def update_fields(self):
        """ MAJ des champs par appel des propriétés du maillage"""
        self.__champs_mailles = \
            {CellPositionField: self.__mesh_instance.coord_elements_field,
             PressureField: self.__mesh_instance.pressure_t_field,
             DensityField: self.__mesh_instance.rho_t_field,
             InternalEnergyField: self.__mesh_instance.nrj_t_field,
             PseudoViscosityField: self.__mesh_instance.pseudo_field
            }

    def create_figure_for_cell_field(self, field_X, field_Y):
        """
        Création des figures pour les champs aux mailles
        (l'axe des X est donc l'abscisse des mailles) 
        """
        try:
            X = np.array(self.__champs_mailles[field_X])
        except ValueError as ve:
            print "Le champ {} est inconnu!".format(field_X)
            raise ve
        try:
            Y = np.array(self.__champs_mailles[field_Y])
        except ValueError as ve:
            print "Le champ {} est inconnu!".format(field_Y)
            raise ve
        if self._dump:
            phyfig = PhysicFigure(X, Y, xlabel=field_X.label, ylabel=field_Y.label,
                              titre=field_Y.titre, save_path=field_Y.results_path)
        else:
            phyfig = PhysicFigure(X, Y, xlabel=field_X.label, ylabel=field_Y.label,
                              titre=field_Y.titre)
        phyfig.set_y_limit(field_Y.val_min, field_Y.val_max)
        phyfig.set_x_limit(field_X.val_min, field_X.val_max)
        return phyfig

    def populate_figs(self):
        """
        Création des figures associées à chacun des champs et ajout à 
        la liste des figures
        """
        champ_X = CellPositionField
        for champ_Y in self.__champs_mailles.keys():
            if champ_Y != champ_X:
                fig = self.create_figure_for_cell_field(champ_X, champ_Y)
                self.__figures_mailles.append((fig, champ_X, champ_Y))
        if self._show:
            plt.show(block=False)
        if self._dump:
            self.create_reps()

    def update_figs(self, title_compl=None):
        """
        MAJ des champs puis des figures correspondantes
        """
        self.update_fields()
        for (fig, champ_x, champ_y) in self.__figures_mailles:
            champ_x_valeurs = self.__champs_mailles[champ_x]
            champ_y_valeurs = self.__champs_mailles[champ_y]
            fig.update(champ_x_valeurs, champ_y_valeurs, title_compl)
        if self._show:
            plt.show(block=False)

    def create_reps(self):
        """
        Création des répertoires où sont stockées les figures
        """
        for (_, _, field) in self.__figures_mailles:
            path = field.results_path
            if (exists(path)):
                msg = "Le chemin {:s} existe déjà!"
                msg += "\nAbandon pour éviter d'écraser des données!"
                raise SystemExit(msg)
            else:
                makedirs(path)