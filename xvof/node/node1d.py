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
    def __init__(self, poz_init=np.zeros(1), vit_init=np.zeros(1),
                section=1.):

        Node.__init__(self, dim=1, position_initiale=poz_init,
                      vitesse_initiale=vit_init)

        #
        self._section = section
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

    def calculer_nouvo_force(self, elements_voisins, isleftboundary=False, isrightboundary=False):
        """
        Calcul de la force agissant sur le noeud
        """
        if isrightboundary:
            pgauche = elements_voisins[0].pressure.new_value + \
            elements_voisins[0].pseudo.current_value
            self._force[:] = pgauche * self.section
        elif isleftboundary:
            pdroite = elements_voisins[0].pressure.new_value + \
            elements_voisins[0].pseudo.current_value
            self._force[:] = -pdroite * self.section
        else:
            pgauche = elements_voisins[0].pressure.new_value + \
                elements_voisins[0].pseudo.current_value
            pdroite = elements_voisins[1].pressure.new_value + \
                elements_voisins[1].pseudo.current_value
            self._force[:] = (pgauche - pdroite) * self.section

    def calculer_nouvo_vitesse(self, delta_t):
        """
        Calcul de la vitesse au demi pas de temps supérieur
        """
#        self._upundemi = self.force / self.masse * delta_t + self.umundemi
        self._upundemi = self.umundemi + self.force * self.invmasse * delta_t

    def appliquer_pression(self, pression):
        """
        Appliquer une pression sur le noeud
        """
        self._force[:] += pression * self.section


if __name__ == "__main__":
    MY_NODE = Node1d(123, section=1.0e-06)
    MY_NODE.infos()
