# -*- coding: iso-8859-1 -*-
"""
Implementing class cell
"""

import numpy as np
from abc import abstractmethod
from copy import deepcopy

from xvof.data.data_container import DataContainer
from xvof.fields.field import Field
from xvof.fields.fieldsmanager import FieldManager


class Cell(object):
    """
    Un objet Element représente l'ensemble des éléments du maillage.
    Ses différents membres sont essentiellement des vecteurs de nbr_of_cells lignes.

    L'organisation en mémoire est comme en C/C++ c'est à dire 'row wise'. C'est pour cette raison
    que les lignes de chacun des vecteurs représentent les mailles. Ce faisant on a par exemple
    toutes les pressions des mailles contigues en mêmoire. (vectorisation, localisation spatiale)
    """

    @classmethod
    def get_coordinates(cls, nbr_cells, topologie, vecteur_x_node, vecteur_y_node=None, vecteur_z_node=None):
        """
        Détermine le vecteur position du centre de l'élément au temps t

        :param nbr_cells: number of cells
        :param topologie: topologie du calcul
        :param vecteur_x_node: vecteur des coordonnées x des noeuds
        :param vecteur_y_node: vecteur des coordonnées y des noeuds
        :param vecteur_z_node: vecteur des coordonnées z des noeuds

        :type topologie: Topology
        :type vecteur_x_node: numpy.array([nbr_of_nodes, 1], dtype=np.float64, order='C')
        :type vecteur_y_node: numpy.array([nbr_of_nodes, 1], dtype=np.float64, order='C')
        :type vecteur_z_node: numpy.array([nbr_of_nodes, 1], dtype=np.float64, order='C')
        :return: les vecteurs des coordonnées du centre de chaque élément au temps t
        :rtype:
        """
        vec_coord = np.zeros([nbr_cells, topologie.dimension])

        for ielem in xrange(nbr_cells):
            nodes_index = topologie.getNodesBelongingToCell(ielem)
            vec_coord[ielem][0] = vecteur_x_node[nodes_index].mean()
            if topologie.dimension == 2:
                vec_coord[ielem][1] = vecteur_y_node[nodes_index].mean()
            if topologie.dimension == 3:
                vec_coord[ielem][2] = vecteur_z_node[nodes_index].mean()
        return vec_coord

    def __init__(self, number_of_elements):
        self._shape = [number_of_elements, ]
        self._dt = np.zeros(self._shape, dtype=np.float64, order='C')
        self._size_t = np.zeros(self._shape, dtype=np.float64, order='C')
        self._size_t_plus_dt = np.zeros(self._shape, dtype=np.float64, order='C')
        self._fields_manager = FieldManager()
        self._fields_manager["Density"] = Field(self._shape[0], DataContainer().material.rho_init,
                                                DataContainer().material.rho_init)
        self._fields_manager["Pressure"] = Field(self._shape[0], DataContainer().material.pression_init,
                                                 DataContainer().material.pression_init)
        self._fields_manager["Energy"] = Field(self._shape[0], DataContainer().material.energie_init,
                                               DataContainer().material.energie_init)
        self._fields_manager["Pseudo"] = Field(self._shape[0])
        self._fields_manager["SoundVelocity"] = Field(self._shape[0])

    @property
    def dt(self):
        """
        Pas de temps critique de la maille
        """
        return self._dt

    @property
    def size_t(self):
        """
        Taille (longueur, aire, volume) de l'élément à l'instant t
        """
        return self._size_t

    @property
    def size_t_plus_dt(self):
        """
        Taille (longueur, aire, volume) de l'élément à l'instant t + dt
        """
        return self._size_t_plus_dt

    @property
    def density(self):
        """
        Champ masse volumique de l'élément
        """
        return self._fields_manager['Density']

    @property
    def pressure(self):
        """
        Champ pression dans l'élément
        """
        return self._fields_manager['Pressure']

    @property
    def sound_velocity(self):
        """
        Champ vitesse du son dans l'élément
        """
        return self._fields_manager['SoundVelocity']

    @property
    def energy(self):
        """
        Champ énergie interne de l'élément
        """
        return self._fields_manager['Energy']

    @property
    def pseudo(self):
        """
        Champ pseudoviscosité dans l'élément
        """
        return self._fields_manager['Pseudo']

    @property
    def fields_manager(self):
        """
        Renvoi une copie du gestionnaire de champs
        """
        return deepcopy(self._fields_manager)

    @property
    def number_of_cells(self):
        return self._shape[0]

    def __str__(self):
        message = "Nombre d'éléments : {:d}".format(self._shape[0])
        return message

    def print_infos(self, index):
        """
        Affichage des informations concernant l'élément d'indice index
        """
        message = "{} {:4d}\n".format(self.__class__, index)
        message += "==> taille à t = {}\n".format(self.size_t[index])
        message += "==> taille à t+dt = {}\n".format(self.size_t_plus_dt[index])
        message += "==> masse volumique à t = {}\n".format(self.density.current_value[index])
        message += "==> masse volumique à t+dt = {}\n". \
            format(self.density.new_value[index])
        message += "==> pression à t = {}\n".format(self.pressure.current_value[index])
        message += "==> pression à t+dt = {}\n". \
            format(self.pressure.new_value[index])
        message += "==> énergie interne à t = {}\n".format(self.energy.current_value[index])
        message += "==> énergie interne à t+dt = {}\n". \
            format(self.energy.new_value[index])
        message += "==> vitesse du son à t = {}\n".format(self.sound_velocity.current_value[index])
        message += "==> vitesse du son à t+dt = {}\n". \
            format(self.sound_velocity.new_value[index])
        print message

    def increment_variables(self):
        """
        Incrémentation des variables
        """
        self._fields_manager.incrementFields()
        self._size_t[:] = self._size_t_plus_dt[:]

    @abstractmethod
    def computeNewPressure(self, mask):
        """
        Algorithme de Newton-Raphson pour déterminer le couple
        energie/pression au pas de temps suivant
        Formulation v-e
        """

    @abstractmethod
    def computeSize(self, topologie):
        """
        Calcul de la taille (longueur, aire, volume) au temps t de l'élément
        """

    @abstractmethod
    def computeNewSize(self, topologie, mask, time_step=None):
        """
        Calcul de la nouvelle taille (longueur, aire, volume) de l'élément
        """

    @abstractmethod
    def computeNewDensity(self, mask):
        """
        Calcul de la densité à l'instant t+dt basé sur
        la conservation de la masse
        """

    @abstractmethod
    def computeNewPseudo(self, time_step, mask):
        """
        Calcul de la nouvelle pseudo
        """

    @abstractmethod
    def computeNewTimeStep(self, mask):
        """
        Calcul du nouveau pas de temps
        """
