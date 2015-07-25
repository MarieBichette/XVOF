#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de base définissant un élément
"""

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# ########### IMPORTATIONS DIVERSES  ####################
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
from abc import abstractmethod
from copy import deepcopy

import numpy as np
from xvof.fields.fieldsmanager import FieldManager


class Element(object):
    """
    Une classe pour les éléments
    """
    @classmethod
    def getCoordinates(cls, noeuds):
        """
        Position du centre de l'élément au temps t
        """
        vec_coord = np.zeros(noeuds[0].dimension)
        for nod in noeuds:
            vec_coord += nod.coordt
        return vec_coord / len(noeuds)

    def __init__(self, proprietes):
        self._index = -1
        self._dt = 0.
        self._size_t = 0.
        self._size_t_plus_dt = 0.
        self._properties = proprietes
        self._fields_manager = FieldManager()
        self._fields_manager.addClassicalField('Density', proprietes.material.rho_init,
                                               proprietes.material.rho_init)
        self._fields_manager.addClassicalField('Pressure', proprietes.material.pression_init,
                                               proprietes.material.pression_init)
        self._fields_manager.addClassicalField('Pseudo')
        self._fields_manager.addClassicalField('SoundVelocity')
        self._fields_manager.addClassicalField('Energy', proprietes.material.energie_init,
                                               proprietes.material.energie_init)
    ##############################################################
    # DEFINITIONS DES PROPRIETES
    ##############################################################
    #

    @property
    def index(self):
        """
        Indice global de l'élément
        """
        return self._index

    @index.setter
    def index(self, index):
        """
        Setter de l'indice global de l'élément
        """
        self._index = index

    @property
    def delta_t(self):
        '''
        Pas de temps critique de la maille
        '''
        return self._dt

    @property
    def taille_t(self):
        """
        Taille (longueur, aire, volume) de l'élément à l'instant t
        """
        return self._size_t

    @property
    def taille_t_plus_dt(self):
        """
        Taille (longueur, aire, volume) de l'élément à l'instant t + dt
        """
        return self._size_t_plus_dt

    @property
    def proprietes(self):
        """
        Proprietes de l'élément
        """
        return self._properties

    @property
    def density(self):
        """
        Champ masse volumique de l'élément
        """
        return self._fields_manager.getField('Density')

    @property
    def pressure(self):
        """
        Champ pression dans l'élément
        """
        return self._fields_manager.getField('Pressure')

    @property
    def sound_velocity(self):
        """
        Champ vitesse du son dans l'élément
        """
        return self._fields_manager.getField('SoundVelocity')

    @property
    def energy(self):
        """
        Champ énergie interne de l'élément
        """
        return self._fields_manager.getField('Energy')

    @property
    def pseudo(self):
        """
        Champ pseudoviscosité dans l'élément
        """
        return self._fields_manager.getField('Pseudo')

    @property
    def fields_manager(self):
        '''
        Renvoi une copie du gestionnaire de champs
        '''
        return deepcopy(self._fields_manager)

    # ------------------------------------------------------------
    # DEFINITIONS DES METHODES
    # ------------------------------------------------------------
    def __str__(self):
        message = "ELEMENT {:4d} ".format(self._index)
        return message

    def printInfos(self):
        """
        Affichage des informations concernant l'élément
        """
        message = "{} {:4d}\n".format(self.__class__, self._index)
        message += "==> taille à t = {}\n".format(self.taille_t)
        message += "==> taille à t+dt = {}\n".format(self.taille_t_plus_dt)
        message += "==> masse volumique à t = {}\n".format(self.density.current_value)
        message += "==> masse volumique à t+dt = {}\n".\
            format(self.density.new_value)
        message += "==> pression à t = {}\n".format(self.pressure.current_value)
        message += "==> pression à t+dt = {}\n".\
            format(self.pressure.new_value)
        message += "==> énergie interne à t = {}\n".format(self.energy.current_value)
        message += "==> énergie interne à t+dt = {}\n".\
            format(self.energy.new_value)
        message += "==> vitesse du son à t = {}\n".format(self.sound_velocity.current_value)
        message += "==> vitesse du son à t+dt = {}\n".\
            format(self.sound_velocity.new_value)
        print message

    def incrementVariables(self):
        """
        Incrémentation des variables
        """
        self._fields_manager.incrementFields()
        self._size_t = self._size_t_plus_dt
    #############################################################
    # DEFINITIONS DES METHODES VIRTUELLES
    #############################################################

    @abstractmethod
    def computeNewPressure(self):
        """
        Algorithme de Newton-Raphson pour déterminer le couple
        energie/pression au pas de temps suivant
        Formulation v-e
        """

    @abstractmethod
    def computeSize(self, nodes):
        """
        Calcul de la taille (longueur, aire, volume) au temps t de l'élément
        """

    @abstractmethod
    def computeNewSize(self, nodes, time_step=None):
        """
        Calcul de la nouvelle taille (longueur, aire, volume) de l'élément
        """

    @abstractmethod
    def computeNewDensity(self):
        """
        Calcul de la densité à l'instant t+dt basé sur
        la conservation de la masse
        """

    @abstractmethod
    def computeNewPseudo(self, time_step):
        """
        Calcul de la nouvelle pseudo
        """

    @abstractmethod
    def computeNewTimeStep(self):
        """
        Calcul du nouveau pas de temps
        """

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# ######          PROGRAMME PRINCIPAL        ###############
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
if __name__ == "__main__":
    import doctest
    TESTRES = doctest.testmod(verbose=0)
    if TESTRES[0] == 0:
        print "TESTS UNITAIRES : OK"
