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
from xvof.miscellaneous import *

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

    def calculer_nouvo_pression(self):
        """
        Algorithme de Newton-Raphson pour déterminer le couple
        energie/pression au pas de temps suivant
        Formulation v-e

        TEST UNITAIRE
        >>> from xvof.miscellaneous import *
        >>> from xvof.equationsofstate import MieGruneisen
        >>> ee = MieGruneisen()
        >>> num_props = numerical_props(0.2, 1.0)
        >>> mat_props = material_props(1.0e+05, 0.0, 8129., ee)
        >>> geom_props = geometrical_props(1.0e-06)
        >>> props = properties(num_props, mat_props, geom_props)
        >>> my_elem = Element(props, 123, 2.5e-03)
        >>> my_elem._rho_t_plus_dt = 9000.0
        >>> my_elem.calculer_nouvo_pression()
        """
        delta_v = 1. / self.rho_t_plus_dt - 1. / self.rho_t
        pression_t = self.pression_t + 2. * self.pseudo
        # Variable du Newton
        nrj_i = self.nrj_t
        # Critère de convergence
        convergence = False
        # Nombre d'itérations
        nit = 0
        #

        def calcul_F_et_dF(enerj):
            """
            Fonction à annuler et sa dérivée pour le schéma VNR
            Formulation v-e
            """
            (p_i, dpsurde, vitson) = \
            self.proprietes.material.eos.solve_ve(1. / self.rho_t_plus_dt, enerj)
            # Fonction à annuler
            func = enerj + p_i * delta_v / 2. + pression_t * delta_v / 2. -\
                self.nrj_t
            # Dérivée de la fonction à annuler
            dfunc = 1 + dpsurde * delta_v / 2.
            return (func, dfunc)
        #
        (fi, dfisurde) = calcul_F_et_dF(nrj_i)
        #
        while(not convergence and (nit < 100)):
            # Correction
            nrj_iplus1 = nrj_i - fi / dfisurde
            nit += 1
            if(abs(fi) < 1e-09):
                convergence = True
                res_nrj = nrj_i
                res_pression_t_plus_dt, dummy, res_cson = \
                self.proprietes.material.eos.solve_ve(
                    1. / self.rho_t_plus_dt, res_nrj)
                break
            # Incrémentation
            nrj_i = nrj_iplus1
            #
            (fi, dfisurde) = calcul_F_et_dF(nrj_i)
            if(abs(dfisurde) < 1.e-09):
                print "Sortie du NR par manque de pente :-)"
                convergence = True
                res_nrj = nrj_i
                res_pression_t_plus_dt = p_i
                res_cson = cel_son
                break
        if(nit == 100):
            print "Erreur de convergence du NR"
            print "fi=", fi
            print "nit=", nit
            exit(255)
        return res_nrj, res_pression_t_plus_dt, res_cson

    #------------------------------------------------------------
    # DEFINITIONS DES METHODES VIRTUELLES
    #------------------------------------------------------------
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

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#######          PROGRAMME PRINCIPAL        ###############
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
if __name__ == "__main__":
    import doctest
    testres = doctest.testmod(verbose=0)
    if(testres[0] == 0):
        print "TESTS UNITAIRES : OK"