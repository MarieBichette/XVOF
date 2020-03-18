#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Module définissant la classe Node1d
"""
from xvof.src.mass_matrix.mass_matrix_utilities import multiplicationMasse
from xvof.src.data.data_container import DataContainer
from xvof.src.node import Node
import numpy as np


class OneDimensionNode(Node):
    """
    Un objet Node1d représente l'ensemble des noeuds 1d du maillage
    """
    def __init__(self, nbr_of_nodes, poz_init, vit_init,
                 section=1.):

        super(OneDimensionNode, self).__init__(nbr_of_nodes, position_initiale=poz_init,
                                                dim=1, vitesse_initiale=vit_init)
        self._section = section
        self._nbr_of_nodes = nbr_of_nodes
        # Par définition, ces noeuds ne sont pas enrichis
        self._classical = np.empty([self.number_of_nodes, ], dtype=np.bool, order='C')
        self._classical[:] = True
        self._enrichment_not_concerned = np.copy(self._classical)

        self.nodes_in_projectile = poz_init[:, 0] <= DataContainer().geometric.initial_interface_position
        self.nodes_in_target = poz_init[:, 0] >= DataContainer().geometric.initial_interface_position
        # on considère que le noeud de contact appartient à la fois au projectile et à la cible

    @property
    def classical(self):
        """
        :return: boolean mask indicating which nodes are classical
        """
        return self._classical

    @property
    def enrichment_not_concerned(self):
        return self._enrichment_not_concerned

    @property
    def velocity_field(self):
        """
        Champ de vitesse vraie
        """
        return self._upundemi

    @property
    def section(self):
        """
        Surface associée au noeud
        """
        return self._section

    def infos(self, index):
        """
        Affichage des informations concernant le noeud d'indice index

        :param index: indice du noeud à afficher
        :type index: int
        """
        Node.infos(self, index)
        message = "==> section = {:5.4g}".format(self.section)
        print message

    def compute_new_force(self, topologie, contrainte):
        """
        Calcul des forces agissant sur les noeuds

        :param topologie: topologie du calcul
        :param contrainte : tenseur des contriante de cauchy sigma xx

        :type topologie: Topology
        :type contrainte: numpy.array([nbr_of_node-1, 1], dtype=np.float64, order='C')
        """
        # Suppose les éléments voisins triés par position croissante
        connectivity = topologie.cells_in_contact_with_node[1:-1]
        self._force[1:-1][:, 0] = (contrainte[connectivity][:, 1] - contrainte[connectivity][:, 0]) * self.section

        ind_node = 0
        elements_voisins = topologie.getCellsInContactWithNode(ind_node)
        self._force[ind_node] = contrainte[elements_voisins][0] * self.section

        ind_node = self.number_of_nodes - 1
        elements_voisins = topologie.getCellsInContactWithNode(ind_node)
        self._force[ind_node] = -contrainte[elements_voisins][0] * self.section

    def compute_new_velocity(self, delta_t, mask, matrice_masse):
        """
        Calcul de la vitesse au demi pas de temps supérieur

        :param delta_t: pas de temps
        :type delta_t: float
        :param mask : noeuds sélectionnés pour calculer avec cette méthode
            (typiquement : OneDimensionEnrichedNode.classical / .enriched )
            s'applique sur les degrés de liberté classiques des noeuds classiques ou enrichis
            ne prend pas en compte le couplage entre degrés de liberté enrichis/classiques
        :type mask : tableau de booléens
        :param matrice_masse : inverse de la matrice de masse
        """
        # ddl classique de noeud classique (sauf 0 1 2 3 quand enrichissement)
        # = noeuds classiques non concernés par l'enrichissement
        self._upundemi[mask] = self._umundemi[mask] + multiplicationMasse(matrice_masse, self.force[mask]) * delta_t

    def compute_complete_velocity_field(self):
        """
        Calcul du champ de vitesse vraie
        """
        self._v_field = np.copy(self._upundemi)

    def apply_correction_reference_bar(self, delta_t, inv_complete_mass_matrix,
                                       inv_wilkins_mass_matrix, mask):
        """
        Apply a correction on velocity field to compute velocity from exact(non lumped) mass matrix for elements in mask
        """
        self._upundemi[mask] += multiplicationMasse(inv_complete_mass_matrix, self.force[mask]) * delta_t
        self._upundemi[mask] -= multiplicationMasse(inv_wilkins_mass_matrix, self.force[mask]) * delta_t

    def apply_pressure(self, ind_node, pressure):
        """
        Appliquer une pressure sur le noeud d'indice "ind_node"

        :param ind_node: indice du noeud sur lequel appliquer la pressure
        :param pressure: pressure à appliquer

        :type ind_node: int
        :type pressure: float
        """
        self._force[ind_node] += pressure * self.section

    def apply_velocity_boundary_coundition(self, ind_node, velocity):
        """
        Appliquer une CL en vitesse sur le noeud ind_node
        :param ind_node:
        :param velocity:
        :return:
        """
        self._upundemi[ind_node] = velocity

