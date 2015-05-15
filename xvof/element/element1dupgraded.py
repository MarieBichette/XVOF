#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant un élément enrichi en 1d
"""
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
############ IMPORTATIONS DIVERSES  ####################
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
import numpy as np
from xvof.element import Element1d
from xvof.node import Node1dUpgraded


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
####### DEFINITION DES CLASSES & FONCTIONS  ###############
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
class Element1dUpgraded(Element1d):
    """
    Une classe pour les éléments enrichis dans le cas 1d
    """
    @classmethod
    def from_geom_to_enrich(cls, champ_gauche, champ_droite):
        """
        Renvoi le champ enrichi à partir des champs gauche et droite
        """
        return (champ_droite - champ_gauche) * 0.5

    @classmethod
    def from_geom_to_classic(cls, champ_gauche, champ_droite):
        """
        Renvoi le champ classique à partir des champs gauche et droite
        """
        return (champ_droite + champ_gauche) * 0.5

    @classmethod
    def from_enrich_to_gauche(cls, champ_classic, champ_enrich):
        """
        Renvoi le champ à gauche d'après les champs classsique et enrichis
        """
        return (champ_classic - champ_enrich)

    @classmethod
    def from_enrich_to_droite(cls, champ_classic, champ_enrich):
        """
        Renvoi le champ à droite d'après les champs classsique et enrichis
        """
        return (champ_classic + champ_enrich)

    def __init__(self, element_origin, pos_discontin):
        Element1d.__init__(self, element_origin.proprietes,
                           element_origin.indice, element_origin.noeuds)
        #
        if(pos_discontin < 0.) or (pos_discontin > 1.):
            message = "La position de la discontinuité dans"
            message += " l'élément enrichi doit être comprise entre 0 et 1!"
            raise SystemExit(message)
        # Les noeuds d'un élément enrichi sont également enrichis
        self._noeuds = map(Node1dUpgraded, self.noeuds)
        #
        self._pression_t = element_origin.pression_t
        self._pression_t_plus_dt = element_origin.pression_t_plus_dt
        self._rho_t = element_origin.rho_t
        self._rho_t_plus_dt = element_origin.rho_t_plus_dt
        self._nrj_t = element_origin.nrj_t
        self._nrj_t_plus_dt = element_origin.nrj_t_plus_dt
        self._pseudo_plus_un_demi = element_origin.pseudo
        self._cson_t = element_origin.cson_t
        self._cson_t_plus_dt = element_origin.cson_t_plus_dt
        #
        self._pression_t_enrichi = 0.
        self._pression_t_plus_dt_enrichi = 0.
        self._rho_t_enrichi = 0.
        self._rho_t_plus_dt_enrichi = 0.
        self._nrj_t_enrichi = 0.
        self._nrj_t_plus_dt_enrichi = 0.
        self._pseudo_plus_un_demi_enrichi = 0.
        self._cson_t_enrichi = 0.
        self._cson_t_plus_dt_enrichi = 0.
        #
        self._taille_gauche_t = element_origin.taille_t * pos_discontin
        self._taille_gauche_t_plus_dt = \
            element_origin.taille_t_plus_dt * pos_discontin
        self._taille_droite_t = element_origin.taille_t * (1. - pos_discontin)
        self._taille_droite_t_plus_dt = \
            element_origin.taille_t_plus_dt * (1. - pos_discontin)
        #
        self._to_gauche = Element1dUpgraded.from_enrich_to_gauche
        self._to_droite = Element1dUpgraded.from_enrich_to_droite
        self._to_classic = Element1dUpgraded.from_geom_to_classic
        self._to_enrich = Element1dUpgraded.from_geom_to_enrich

    @property
    def taille_t_gauche(self):
        """
        Taille de la partie gauche de l'élément au temps t
        """
        return self._taille_gauche_t

    @property
    def taille_t_droite(self):
        """
        Taille de la partie droite de l'élément au temps t
        """
        return self._taille_droite_t

    @property
    def taille_t_plus_dt_gauche(self):
        """
        Taille de la partie gauche de l'élément au temps t+dt
        """
        return self._taille_gauche_t_plus_dt

    @property
    def taille_t_plus_dt_droite(self):
        """
        Taille de la partie droite de l'élément au temps t+dt
        """
        return self._taille_droite_t_plus_dt

    @property
    def pression_t_gauche(self):
        """
        Pression dans la partie gauche de l'élément au temps t
        """
        return self._to_gauche(self._pression_t,
            self._pression_t_enrichi)

    @property
    def pression_t_droite(self):
        """
        Pression dans la partie droite de l'élément au temps t
        """
        return self._to_droite(self._pression_t,
            self._pression_t_enrichi)

    @property
    def rho_t_gauche(self):
        """
        Densité dans la partie gauche de l'élément au temps t
        """
        return self._to_gauche(self._rho_t, self._rho_t_enrichi)

    @property
    def rho_t_droite(self):
        """
        Densité dans la partie droite de l'élément au temps t
        """
        return self._to_droite(self._rho_t, self._rho_t_enrichi)

    @property
    def rho_t_plus_dt_gauche(self):
        """
        Densité dans la partie gauche de l'élément au temps t+dt
        """
        return self._to_gauche(self._rho_t_plus_dt,
                               self._rho_t_plus_dt_enrichi)

    @property
    def rho_t_plus_dt_droite(self):
        """
        Densité dans la partie droite de l'élément au temps t+dt
        """
        return self._to_droite(self._rho_t_plus_dt,
                               self._rho_t_plus_dt_enrichi)

    @property
    def nrj_t_gauche(self):
        """
        Densité dans la partie gauche de l'élément au temps t
        """
        return self._to_gauche(self._nrj_t, self._nrj_t_enrichi)

    @property
    def nrj_t_droite(self):
        """
        Energie dans la partie droite de l'élément au temps t
        """
        return self._to_droite(self._nrj_t, self._nrj_t_enrichi)

    @property
    def nrj_t_plus_dt_gauche(self):
        """
        Energie dans la partie gauche de l'élément au temps t+dt
        """
        return self._to_gauche(self._nrj_t_plus_dt,
                               self._nrj_t_plus_dt_enrichi)

    @property
    def nrj_t_plus_dt_droite(self):
        """
        Energie dans la partie droite de l'élément au temps t+dt
        """
        return self._to_droite(self._nrj_t_plus_dt,
                               self._nrj_t_plus_dt_enrichi)

    @property
    def cson_t_gauche(self):
        """
        Vitesse du son dans la partie gauche de l'élément au temps t
        """
        return self._to_gauche(self._cson_t, self._cson_t_enrichi)

    @property
    def cson_t_droite(self):
        """
        Vitesse du son dans la partie droite de l'élément au temps t
        """
        return self._to_droite(self._cson_t, self._cson_t_enrichi)

    @property
    def cson_t_plus_dt_gauche(self):
        """
        Vitesse du son dans la partie gauche de l'élément au temps t+dt
        """
        return self._to_gauche(self._cson_t_plus_dt,
                               self._cson_t_plus_dt_enrichi)

    @property
    def cson_t_plus_dt_droite(self):
        """
        Vitesse du son dans la partie droite de l'élément au temps t+dt
        """
        return self._to_droite(self._cson_t_plus_dt,
                               self._cson_t_plus_dt_enrichi)

    @property
    def pseudo_gauche(self):
        """
        Pseudo viscosité dans la partie gauche de l'élément
        """
        return self._to_gauche(self._pseudo_plus_un_demi,
                               self._pseudo_plus_un_demi_enrichi)

    @property
    def pseudo_droite(self):
        """
        Pseudo viscosité dans la partie droite de l'élément
        """
        return self._to_droite(self._pseudo_plus_un_demi,
                               self._pseudo_plus_un_demi_enrichi)

    # --------------------------------------------------------
    #        DEFINITION DES METHODES                         #
    # --------------------------------------------------------
    def __str__(self):
        message = "ELEMENT ENRICHI {:4d} ".format(self.indice)
        return message

    def infos(self):
        """
        Affichage des informations concernant l'élément
        """
        Element1d.infos(self)
        message = "==> masse volumique à gauche à t = {}\n".\
            format(self.rho_t_gauche)
        message += "==> masse volumique à droite à t = {}\n".\
            format(self.rho_t_droite)
        message += "==> masse volumique à gauche à t+dt = {}\n".\
            format(self.rho_t_plus_dt_gauche)
        message += "==> masse volumique à droite à t+dt = {}\n".\
            format(self.rho_t_plus_dt_droite)
        message += "==> taille à gauche à t = {}\n".\
            format(self.taille_t_gauche)
        message += "==> taille à droite à t = {}\n".\
            format(self.taille_t_droite)
        message += "==> taille à gauche à t+dt = {}\n".\
            format(self.taille_t_plus_dt_gauche)
        message += "==> taille à droite à t+dt = {}\n".\
            format(self.taille_t_plus_dt_droite)
        message += "==> pression à gauche à t = {}\n".\
            format(self.pression_t_gauche)
        message += "==> pression à droite à t = {}\n".\
            format(self.pression_t_droite)
        message += "==> pression à gauche à t+dt = {}\n".\
            format(self.pression_t_plus_dt_gauche)
        message += "==> pression à droite à t+dt = {}\n".\
            format(self.pression_t_plus_dt_droite)
        message += "==> vitesse du son à gauche à t = {}\n".\
            format(self.cson_t_gauche)
        message += "==> vitesse du son à droite à t = {}\n".\
            format(self.cson_t_droite)
        message += "==> vitesse du son à gauche à t+dt = {}\n".\
            format(self.cson_t_plus_dt_gauche)
        message += "==> vitesse du son à droite à t+dt = {}\n".\
            format(self.cson_t_plus_dt_droite)
        message += "==> énergie à gauche à t = {}\n".\
            format(self.nrj_t_gauche)
        message += "==> énergie à droite à t = {}\n".\
            format(self.nrj_t_droite)
        message += "==> énergie à gauche à t+dt = {}\n".\
            format(self.nrj_t_plus_dt_gauche)
        message += "==> énergie à droite à t+dt = {}\n".\
            format(self.nrj_t_plus_dt_droite)
        message += "==> pseudo à gauche = {}\n".\
            format(self.pseudo_gauche)
        message += "==> pseudo à droite = {}\n".\
            format(self.pseudo_droite)
        print message

    def calculer_nouvo_pression(self):
        """
        Calcul du triplet energie, pression, vitesse du son
        au pas de temps suivant
        Formulation v-e
        """
        nrj_t_plus_dt_g, pression_t_plus_dt_g, cson_t_plus_dt_g = \
            Element1d.newton_raphson_for_ve(self.proprietes.material.eos,
                                            self.rho_t_gauche,
                                            self.rho_t_plus_dt_gauche,
                                            self.pression_t_gauche,
                                            self.pseudo_gauche,
                                            self.nrj_t_gauche)
        nrj_t_plus_dt_d, pression_t_plus_dt_d, cson_t_plus_dt_d = \
            Element1d.newton_raphson_for_ve(self.proprietes.material.eos,
                                            self.rho_t_droite,
                                            self.rho_t_plus_dt_droite,
                                            self.pression_t_droite,
                                            self.pseudo_droite,
                                            self.nrj_t_droite)
        #
        self._pression_t_plus_dt = \
            self._to_classic(pression_t_plus_dt_g, pression_t_plus_dt_d)
        self._pression_t_plus_dt_enrichi = \
            self._to_enrich(pression_t_plus_dt_g, pression_t_plus_dt_d)
        #
        self._nrj_t_plus_dt = \
            self._to_classic(nrj_t_plus_dt_g, nrj_t_plus_dt_d)
        self._nrj_t_plus_dt_enrichi = \
            self._to_enrich(nrj_t_plus_dt_g, nrj_t_plus_dt_d)
        #
        self._cson_t_plus_dt = \
            self._to_classic(cson_t_plus_dt_g, cson_t_plus_dt_d)
        self._cson_t_plus_dt_enrichi = \
            self._to_enrich(cson_t_plus_dt_g, cson_t_plus_dt_d)

    def calculer_nouvo_taille(self, delta_t):
        """
        Calcul des nouvelles longueurs de l'élément
        """
        # Les noeuds sont classés par coord croissante
        nod_g = self.noeuds[0]
        nod_d = self.noeuds[1]
        self._taille_gauche_t_plus_dt = self.taille_t_gauche + \
            (0.5 * (nod_d.upundemi_classique - nod_g.upundemi_enrichi) - 
             0.5 * (nod_g.upundemi_classique - nod_g.upundemi_enrichi)) \
            * delta_t
        self._taille_droite_t_plus_dt = self.taille_t_droite + \
            (0.5 * (nod_d.upundemi_classique + nod_d.upundemi_enrichi) - 
             0.5 * (nod_g.upundemi_classique + nod_d.upundemi_enrichi)) \
            * delta_t

    def calculer_nouvo_densite(self):
        """
        Calcul des nouvelles densités
        """
        densite_gauche_t_plus_dt = self.rho_t_gauche * self.taille_t_gauche \
            / self.taille_t_plus_dt_gauche
        densite_droite_t_plus_dt = self.rho_t_droite * self.taille_t_droite \
            / self.taille_t_plus_dt_droite
        self._rho_t_plus_dt = \
            self._to_classic(densite_gauche_t_plus_dt, densite_droite_t_plus_dt)
        self._rho_t_plus_dt_enrichi = \
            self._to_enrich(densite_gauche_t_plus_dt, densite_droite_t_plus_dt)

    def calculer_nouvo_pseudo(self, delta_t):
        """
        Calcul de la nouvelle pseudo
        """
        pseudo_gauche = \
            Element1d.calculer_pseudo(delta_t, self.rho_t_gauche,
                                      self.rho_t_plus_dt_gauche,
                                      self.taille_t_plus_dt_gauche,
                                      self.cson_t_gauche,
                                      self.proprietes.numeric.a_pseudo,
                                      self.proprietes.numeric.b_pseudo)

        pseudo_droite = \
            Element1d.calculer_pseudo(delta_t, self.rho_t_droite,
                                      self.rho_t_plus_dt_droite,
                                      self.taille_t_plus_dt_droite,
                                      self.cson_t_droite,
                                      self.proprietes.numeric.a_pseudo,
                                      self.proprietes.numeric.b_pseudo)

        self._pseudo_plus_un_demi = \
            self._to_classic(pseudo_gauche, pseudo_droite)
        self._pseudo_plus_un_demi_enrichi = \
            self._to_enrich(pseudo_gauche, pseudo_droite)

    def calculer_nouvo_dt(self):
        """
        Calcul du pas de temps
        """
        cfl = self.proprietes.numeric.cfl
        dt_g = \
            Element1d.calculer_dt(cfl, self.rho_t_gauche,
                                  self.rho_t_plus_dt_gauche,
                                  self.taille_t_plus_dt_gauche,
                                  self.cson_t_plus_dt_gauche,
                                  self.pseudo_gauche)

        dt_d = \
            Element1d.calculer_dt(cfl, self.rho_t_droite,
                                  self.rho_t_plus_dt_droite,
                                  self.taille_t_plus_dt_droite,
                                  self.cson_t_plus_dt_droite,
                                  self.pseudo_droite)

        self._dt = dt_g + dt_d  # Bizarre --> A vérifier

    def incrementer(self):
        """
        Incrémentation des variables
        """
        Element1d.incrementer(self)
        self._pression_t_enrichi = self._pression_t_plus_dt_enrichi
        self._rho_t_enrichi = self._rho_t_plus_dt_enrichi
        self._cson_t_enrichi = self._cson_t_plus_dt_enrichi
        self._nrj_t_enrichi = self._nrj_t_plus_dt_enrichi
