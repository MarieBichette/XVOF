#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de base définissant un élément
"""

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
############ IMPORTATIONS DIVERSES  ####################
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
from abc import abstractmethod
from xvof.miscellaneous import *
import numpy as np
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
####### DEFINITION DES CLASSES & FONCTIONS  ###############
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


class Element(object):
    """
    Une classe pour les éléments
    """
    # pylint: disable-msg=R0902
    # 13 attributs : cela semble raisonnable pour ce cas
    def __init__(self, proprietes, indice, noeuds):
        self._index = indice
        self._dt = 0.
        self._size_t = 0.
        self._size_t_plus_dt = 0.
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
        self.noeuds = noeuds

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
    def delta_t(self):
        """
        Pas de temps de l'élément
        """
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

    @property
    def coord(self):
        """
        Position du centre de l'élément au temps t

        TEST UNITAIRE
        >>> from xvof.node import Node1d
        >>> from xvof.miscellaneous import *
        >>> from xvof.equationsofstate import MieGruneisen
        >>> ee = MieGruneisen()
        >>> num_props = numerical_props(0.2, 1.0, 0.35)
        >>> mat_props = material_props(1.0e+05, 0.0, 8129., ee)
        >>> geom_props = geometrical_props(1.0e-06)
        >>> props = properties(num_props, mat_props, geom_props)
        >>> noda = Node1d(1, np.array([-0.5]))
        >>> nodb = Node1d(2, np.array([0.5]))
        >>> my_elem = Element(props, 1, [noda, nodb])
        >>> my_elem.coord
        array([ 0.])
        """
        vec_coord = np.zeros(self.noeuds[0].dimension)
        for nod in self.noeuds:
            vec_coord += nod.coordt
        return vec_coord / len(self.noeuds)

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
        message = "{} {:4d}\n".format(self.__class__, self.indice)
        message += "==> noeuds = {}\n".format(self.noeuds)
        message += "==> taille à t = {}\n".format(self.taille_t)
        message += "==> taille à t+dt = {}\n".format(self.taille_t_plus_dt)
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

    def incrementer(self):
        """
        Incrémentation des variables
        """
        self._pression_t = self._pression_t_plus_dt
        self._rho_t = self._rho_t_plus_dt
        self._cson_t = self._cson_t_plus_dt
        self._nrj_t = self._nrj_t_plus_dt
        self._size_t_plus_dt = self._size_t
    #------------------------------------------------------------
    # DEFINITIONS DES METHODES VIRTUELLES
    #------------------------------------------------------------
    @abstractmethod
    def calculer_nouvo_pression(self):
        """
        Algorithme de Newton-Raphson pour déterminer le couple
        energie/pression au pas de temps suivant
        Formulation v-e
        """

    @abstractmethod
    def calculer_nouvo_taille(self):
        """
        Calcul de la nouvelle taille (longueur, aire, volume) de l'élément
        """

    @abstractmethod
    def calculer_nouvo_densite(self):
        """
        Calcul de la densité à l'instant t+dt basé sur
        la conservation de la masse
        """

    @abstractmethod
    def calculer_nouvo_pseudo(self):
        """
        Calcul de la nouvelle pseudo
        """

    @abstractmethod
    def calculer_nouvo_dt(self):
        """
        Calcul du nouveau pas de temps
        """

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#######          PROGRAMME PRINCIPAL        ###############
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
if __name__ == "__main__":
    import doctest
    TESTRES = doctest.testmod(verbose=0)
    if(TESTRES[0] == 0):
        print "TESTS UNITAIRES : OK"
