#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant un noeud en 1d
"""
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
############ IMPORTATIONS DIVERSES  ####################
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
import numpy as np
from xvof.node import Node


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
####### DEFINITION DES CLASSES & FONCTIONS  ###############
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
class Node1d(Node):
    """
    Une classe pour les noeuds classiques dans le cas 1d
    """
    def __init__(self, indice, poz_init=np.zeros(1), vit_init=np.zeros(1),
                section=1.):

        Node.__init__(self, dim=1, index=indice, position_initiale=poz_init,
                      vitesse_initiale=vit_init)

        #
        self._section = section
        #
        self._elements_voisins = []
        #
        self._upundemi = vit_init

    #------------------------------------------------------------
    # DEFINITIONS DES PROPRIETES
    #------------------------------------------------------------

    @property
    def section(self):
        """
        Surface associée au noeud
        """
        return self._section

    @property
    def elements_voisins(self):
        """
        Liste des éléments voisins du noeud
        """
        return self._elements_voisins

    @elements_voisins.setter
    def elements_voisins(self, elems):
        """
        Setter des elements voisins. Surcharge de la méthode de Node pour
        s'assurer qu'il n y ait que deux voisins possibles et pour les trier
        de gauche à droite
        """
        self._elements_voisins.extend(elems)
        if(len(self.elements_voisins) > 2):
            message = "En 1d au plus deux éléments peuveut être"
            message += " voisins du {}".format(self)
            message += "\n Liste des elements : {}".format(self.elements_voisins)
            raise SystemExit(message)
        self._elements_voisins = \
        sorted(self._elements_voisins, key=lambda m: m.coord[0])

    #------------------------------------------------------------
    # DEFINITIONS DES METHODES
    #------------------------------------------------------------

    def infos(self):
        """
        Affichage des informations
        """
        Node.infos(self)
        message = "==> section = {:5.4g}".format(self.section)
        print message

    def calculer_nouvo_force(self):
        """
        Calcul de la force agissant sur le noeud
        """
        if(len(self.elements_voisins) == 2):
            pgauche = self.elements_voisins[0].pression_t_plus_dt + \
                self.elements_voisins[0].pseudo
            pdroite = self.elements_voisins[1].pression_t_plus_dt + \
                self.elements_voisins[1].pseudo
            self._force[:] = (pgauche - pdroite) * self.section
        else:
            # Cas des noeuds de bord
            if(self.coordt[0] > self.elements_voisins[0].coord[0]):
                # Noeud du bord droit
                pgauche = self.elements_voisins[0].pression_t_plus_dt + \
                self.elements_voisins[0].pseudo
                self._force[:] = pgauche * self.section
            elif(self.coordt[0] < self.elements_voisins[0].coord[0]):
                # Noeud du bord gauche
                pdroite = self.elements_voisins[0].pression_t_plus_dt + \
                self.elements_voisins[0].pseudo
                self._force[:] = -pdroite * self.section

    def calculer_nouvo_vitesse(self, delta_t):
        """
        Calcul de la vitesse au demi pas de temps supérieur
        """
        self._upundemi = self.force / self.masse * delta_t + self.umundemi

    def appliquer_pression(self, pression):
        """
        Appliquer une pression sur le noeud
        """
        self._force[:] += pression * self.section


if __name__ == "__main__":
    MY_NODE = Node1d(123, section=1.0e-06)
    MY_NODE.infos()
