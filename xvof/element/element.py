#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de base définissant un élément
"""

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
############ IMPORTATIONS DIVERSES  ####################
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
from abc import abstractmethod
import numpy as np
from miscellaneous import *

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
####### DEFINITION DES CLASSES & FONCTIONS  ###############
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


class Element(object):
    """
    Une classe pour les éléments
    """
    def __init__(self, proprietes, indice, taille):
        self._index = indice
        self._size = taille
        self._properties = proprietes
        self._rho_t = proprietes.material.rho_init
        self._rho_t_plus_dt = proprietes.material.rho_init
        self._pression_t = proprietes.material.pression_init
        self._pression_t_plus_dt = proprietes.material.pression_init
        self._pseudo_plus_un_demi = 0.
        self._cson_t = 0.
        self._cson_t_plus_dt = 0.
        self._nrj_t = proprietes.material.energie_init
        self._nrj_t_plus_dt = proprietes.material.energie_init
        self._noeuds = []

    #------------------------------------------------------------
    # DEFINITIONS DES PROPRIETES
    #------------------------------------------------------------
    #
    # Seules les modifications de _noeuds_voisins sont permises
    # Les autres attributs sont accessibles en lecture seule
    #
    @property
    def indice(self):
        """
        Indice global de l'élément
        """
        return self._index

    @property
    def taille(self):
        """
        Taille (longueur, aire, volume) de l'élément
        """
        return self._size

    @property
    def proprietes(self):
        """
        Proprietes de l'élément
        """
        return self._properties

    @property
    def rho_t(self):
        """
        Masse volumique de l'élément au temps t
        """
        return self._rho_t

    @property
    def rho_t_plus_dt(self):
        """
        Masse volumique de l'élément au temps t + dt
        """
        return self._rho_t_plus_dt

    @property
    def pression_t(self):
        """
        Pression dans l'élément au temps t
        """
        return self._pression_t

    @property
    def pression_t_plus_dt(self):
        """
        Pression dans l'élément au temps t + dt
        """
        return self._pression_t_plus_dt

    @property
    def cson_t(self):
        """
        Vitesse du son dans l'élément au temps t
        """
        return self._cson_t

    @property
    def cson_t_plus_dt(self):
        """
        Vitesse du son dans l'élément au temps t + dt
        """
        return self._cson_t_plus_dt

    @property
    def nrj_t(self):
        """
        Energie interne de l'élément au temps t
        """
        return self._nrj_t

    @property
    def nrj_t_plus_dt(self):
        """
        Energie interne dans l'élément au temps t + dt
        """
        return self._nrj_t_plus_dt

    @property
    def pseudo(self):
        """
        Pseudo viscosité dans l'élément
        """
        return self._pseudo_plus_un_demi

    @property
    def noeuds(self):
        """
        Liste des noeuds de l'élément
        """
        return self._noeuds

    @noeuds.setter
    def noeuds(self, node_list):
        """
        Setter des noeuds de l'élément
        """
        self._noeuds[:] = node_list[:]

    #------------------------------------------------------------
    # DEFINITIONS DES METHODES
    #------------------------------------------------------------
    def __str__(self):
        message = "ELEMENT {:4d} ".format(self.indice)
        return message

    def infos(self):
        """
        Affichage des informations concernant l'élément
        """
        message = "{} {:4d}\n".format(self.__class__, self.index)
        message += "==> noeuds = {}\n".format(self.noeuds)
        message += "==> taille = {}\n".format(self.taille)
        message += "==> masse volumique à t = {}\n".format(self.rho_t)
        message += "==> masse volumique à t+dt = {}\n".\
            format(self.rho_t_plus_dt)
        message += "==> pression à t = {}\n".format(self.pression_t)
        message += "==> pression à t+dt = {}\n".\
            format(self.pression_t_plus_dt)
        message += "==> énergie interne à t = {}\n".format(self.nrj_t)
        message += "==> énergie interne à t+dt = {}\n".\
            format(self.nrj_t_plus_dt)
        message += "==> vitesse du son à t = {}\n".format(self.cson_t)
        message += "==> vitesse du son à t+dt = {}\n".\
            format(self.cson_t_plus_dt)
        print message