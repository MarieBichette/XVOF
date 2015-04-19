#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de base définissant un maillage 1d
"""
import numpy as np
from xvof.element.element1d import Element1d
from xvof.equationsofstate.miegruneisen import MieGruneisen
from xvof.miscellaneous import geometrical_props, material_props
from xvof.miscellaneous import numerical_props, properties
from xvof.node.node1d import Node1d


class Mesh1d(object):
    """
    Une classe définissant un maillage 1d
    """
    def __init__(self, proprietes, initial_coordinates=np.linspace(0, 1, 11),
                 initial_velocities=np.zeros(11)):
        if np.shape(initial_coordinates) != np.shape(initial_velocities):
            message = "Les vecteurs initiaux de vitesse et coordonnées"
            message += " n'ont pas la même taille"
            raise ValueError(message)
        if len(np.shape(initial_coordinates)) != 1:
            message = "Il s'agit d'un maillage 1D à initialiser avec des"
            message += " vecteurs à une dimension!"
            raise ValueError(message)
        self.__nbr_nodes = np.shape(initial_coordinates)[0]
        self.__nbr_cells = self.__nbr_nodes - 1
        self.__nodes = []
        self.__cells = []
        # Création des noeuds
        for n in xrange(self.__nbr_nodes):
            poz = initial_coordinates[n]
            vit = initial_velocities[n]
            nod = Node1d(n, poz_init=np.array([poz]),
                         vit_init=np.array([vit]))
            self.__nodes.append(nod)
        # Création des éléments
        for m in xrange(self.__nbr_cells):
            elem = Element1d(proprietes, m, [self.__nodes[m], self.__nodes[m + 1]])
            self.__cells.append(elem)
            self.__nodes[m].elements_voisins = [elem]
            self.__nodes[m + 1].elements_voisins = [elem]
    
    @property
    def nodes(self):
        """ Liste des noeuds """
        return self.__nodes

    @property
    def cells(self):
        """ Liste des éléments """
        return self.__cells

    def calculer_masse_des_noeuds(self):
        """ Calcul de la masse de chaque noeud"""
        for noeud in self.nodes:
            noeud.calculer_masse_wilkins()

    def calculer_nouvo_vit_noeuds(self, delta_t):
        """ Calcul de la nouvelle vitesse de chaque noeud à t+dt"""
        for noeud in self.nodes:
            noeud.calculer_nouvo_vitesse(delta_t)

    def calculer_nouvo_coord_noeuds(self, delta_t):
        """ Calcul des nouvelles coordonnées de chaque noeud à t+dt"""
        for noeud in self.nodes:
            noeud.calculer_nouvo_coord(delta_t)

    def calculer_nouvo_taille_des_elements(self):
        """ Calcul de la nouvelle taille de chaque élément à t+dt"""
        for cell in self.cells:
            cell.calculer_nouvo_taille()

    def calculer_nouvo_densite_des_elements(self):
        """ Calcul des nouvelles densités de chaque élément à t+dt"""
        for cell in self.cells:
            cell.calculer_nouvo_densite()

    def calculer_nouvo_pression_des_elements(self):
        """ Calcul des nouvelles pressions de chaque élément à t+dt"""
        for cell in self.cells:
            cell.calculer_nouvo_pression()

    def calculer_nouvo_force_des_noeuds(self):
        """ Calcul des nouvelles forces de chaque noeud à t+dt"""
        for noeud in self.nodes:
            noeud.calculer_nouvo_force()

    def incrementer(self):
        """ Passage au pas de temps suivant"""
        for noeud in self.nodes:
            noeud.incrementer()
        for cell in self.cells:
            cell.incrementer()