#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Module définissant la classe Node1d
"""
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# ########### IMPORTATIONS DIVERSES  ####################
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
import numpy as np
from xvof.node import Node


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# ###### DEFINITION DES CLASSES & FONCTIONS  ###############
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
class Node1d(Node):
    """
    Un objet Node1d représente l'ensemble des noeuds 1d du maillage
    """
    def __init__(self, nbr_of_nodes, poz_init, vit_init,
                 section=1.):

        Node.__init__(self, nbr_of_nodes, position_initiale=poz_init,
		      dim=1, vitesse_initiale=vit_init)
        self._section = section

    # ------------------------------------------------------------
    # DEFINITIONS DES PROPRIETES
    # ------------------------------------------------------------

    @property
    def section(self):
        """
        Surface associée au noeud
        """
        return self._section

    # ------------------------------------------------------------
    # DEFINITIONS DES METHODES
    # ------------------------------------------------------------

    def infos(self, index):
        """
        Affichage des informations concernant le noeud d'indice index

        :param index: indice du noeud à afficher
        :type index: int
        """
        Node.infos(self, index)
        message = "==> section = {:5.4g}".format(self.section)
        print message

    def calculer_nouvo_force(self, topologie, vecteur_pression_maille, vecteur_pseudo_maille):
        """
        Calcul des forces agissant sur les noeuds

        :param topologie: topologie du calcul
        :param vecteur_pression_maille: vecteur des pressions de chaque élément + 2 pressions nulles à gauche et à droite
        :param vecteur_pseudo_maille: vecteur des pseudoviscosite de chaque élément

        :type topologie: Topology
        :type vecteur_pression_maille: numpy.array([nbr_of_nodes+2, 1], dtype=np.float64, order='C')
        :type vecteur_pseudo_maille: numpy.array([nbr_of_nodes, 1], dtype=np.int64, order='C')
        """
        # Suppose les éléments voisins triés par position croissante
        connectivity = np.array(topologie._cells_in_contact_with_node[1 : self.number_of_nodes - 1])
        pressure = self.section * (vecteur_pression_maille[connectivity] + vecteur_pseudo_maille[connectivity])
        self._force[1:self.number_of_nodes - 1][:, 0] = (pressure[:, 0] - pressure[:, 1]).reshape(self.number_of_nodes - 2)
        ind_node = 0
        elements_voisins = topologie.getCellsInContactWithNode(ind_node)
        pressure = self.section * (vecteur_pression_maille[elements_voisins] + vecteur_pseudo_maille[elements_voisins])
        self._force[ind_node] = -pressure[0]
        ind_node = self.number_of_nodes - 1
        elements_voisins = topologie.getCellsInContactWithNode(ind_node)
        pressure = self.section * (vecteur_pression_maille[elements_voisins] + vecteur_pseudo_maille[elements_voisins])
        self._force[ind_node] = pressure[0]

    def calculer_nouvo_vitesse(self, delta_t):
        """
        Calcul de la vitesse au demi pas de temps supérieur

        :param delta_t: pas de temps
        :type delta_t: float
        """
        self._upundemi = self.umundemi + self.force * self.invmasse * delta_t

    def applyPressure(self, ind_node, pressure):
        """
        Appliquer une pressure sur le noeud d'indice "ind_node"

        :param ind_node: indice du noeud sur lequel appliquer la pressure
        :param pressure: pressure à appliquer

        :type ind_node: int
        :type pressure: float
        """
        self._force[ind_node] += pressure * self.section

