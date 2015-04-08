#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant un noeud enrichi en 1d
"""
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
############ IMPORTATIONS DIVERSES  ####################
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
from xvof.node import Node1d
import numpy as np

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
####### DEFINITION DES CLASSES & FONCTIONS  ###############
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


class Node1dUpgraded(Node1d):
    """
    Une classe pour les noeuds enrichis dans le cas 1d
    """
    def __init__(self, indice, poz_init=np.zeros(1), vit_init=np.zeros(1),
                section=1.):
        Node1d.__init__(self, dim=1, index=indice, position_initiale=poz_init,
                      vitesse_initiale=vit_init, section=section)

        self._umundemi_classique = vit_init[:]
        self._upundemi_classique = vit_init[:]
        self._force_classique = np.zeros(1, dtype=float)
        #==> Toutes les variables enrichies sont initialisées à 0
        self._umundemi_enrichi = np.zeros(1, dtype=float)
        self._upundemi_enrichi = np.zeros(1, dtype=float)
        self._force_enrichi = np.zeros(1, dtype=float)
        #
        self.__position_relative = None

    #------------------------------------------------------------
    # DEFINITIONS DES PROPRIETES
    #------------------------------------------------------------

    @property
    def position_relative(self):
        """
        Position du noeud par rapport à la discontinuité
        """
        return self.__position_relative

    @position_relative.setter
    def position_relative(self, pos):
        """
        Setter de la position relative
        """
        if (pos not in (-1, 1)):
            message = "La position relative du noeud ne peut être que :\n"
            message += " -1 si il est à gauche de la discontinuité\n"
            message += " +1 si il est à droite"
            raise SystemExit(message)
        self.__position_relative = pos

    @property
    def umundemi_classique(self):
        return self._umundemi_classique

    @property
    def umundemi_enrichi(self):
        return self._umundemi_enrichi

    @property
    def upundemi_classique(self):
        return self._upundemi_classique

    @property
    def upundemi_enrichi(self):
        return self._upundemi_enrichi

    @property
    def force_classique(self):
        return self._force_classique

    @property
    def force_enrichi(self):
        return self._force_enrichi

    #------------------------------------------------------------
    # DEFINITIONS DES METHODES
    #------------------------------------------------------------

    def infos(self):
        """
        Affichage des informations
        """
        Node1d.infos(self)
        message = "==> vitesse classique à t-1/2 = {}\n".\
        format(self.umundemi_classique)
        message += "==> vitesse enrichie à t-1/2 = {}\n".\
        format(self.umundemi_enrichi)
        message += "==> vitesse classique à t+1/2 = {}\n".\
        format(self.upundemi_classique)
        message += "==> vitesse enrichie à t+1/2 = {}\n".\
        format(self.upundemi_enrichi)
        message += "==> force classique à t-1/2 = {}\n".\
        format(self.force_classique)
        message += "==> force enrichie à t-1/2 = {}\n".\
        format(self.force_enrichi)
        message += "==> position relative  = {:2d}".\
        format(self.position_relative)

    def initialize(self, vitesse_t_m, vitesse_t_p, force):
        """
        Initialisation des champs classiques

        @param vitesse_t_m : vecteur vitesse au demi pas de temps précédent
        @param vitesse_t_p : vecteur vitesse au demi pas de temps suivant
        @param force : vecteur force
        """
        self._umundemi_classique = vitesse_t_m[:]
        self._upundemi_classique = vitesse_t_p[:]
        self._force_classique = force[:]

    def calculer_nouvo_vitesse(self, delta_t):
        """
        Calcul de la vitesse au demi pas de temps supérieur
        """
        self._upundemi_enrichi =\
            self.force_enrichi / self.masse * delta_t + self.umundemi_enrichi
        self._upundemi_classique =\
            self.force_classique / self.masse * delta_t + \
                self.umundemi_classique
        self._upundemi =\
            self.upundemi_classique + \
            self.position_relative * self.upundemi_enrichi

    def calculer_nouvo_force(self):
        """
        Calcul de la force agissant sur le noeud
        """
        self._force_classique[:] = (self.elements_voisins[0].pressure -
            self.elements_voisins[1].pressure) * self.section
        self._force_enrichi[:] = (-self.elements_voisins[0].pressure -
            self.elements_voisins[1].pressure) * self.section
        self._force = None

    def incrementer(self):
        """
        Mise à jour de la vitesse et de la coordonnée du noeud
        pour passer au pas de temps suivant.
        """
        Node1d.incrementer(self)
        self._umundemi_classique[:] = self.upundemi_classique[:]
        self._umundemi_enrichi[:] = self.upundemi_enrichi[:]