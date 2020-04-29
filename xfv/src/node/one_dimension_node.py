#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Module définissant la classe Node1d
"""
import numpy as np
from xfv.src.mass_matrix.mass_matrix_utilities import multiplication_masse
from xfv.src.data.data_container import DataContainer
from xfv.src.node import Node


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
        # By definition, these nodes are not enriched
        self._classical = np.empty([self.number_of_nodes, ], dtype=np.bool, order='C')
        self._classical[:] = True
        self._enrichment_not_concerned = np.copy(self._classical)

        interface_position = DataContainer().geometric.initial_interface_position
        self.nodes_in_projectile = poz_init[:, 0] <= interface_position
        self.nodes_in_target = poz_init[:, 0] >= interface_position
        # the node in contact belongs both to the target and projectile

        # Velocity field that takes into account the enriched degrees of freedom influence
        self._v_field = np.zeros_like(self.upundemi)

    @property
    def classical(self):
        """
        :return: boolean mask indicating which nodes are classical
        """
        return self._classical

    @property
    def enrichment_not_concerned(self):
        """
        :return: boolean mask indicating which nodes are concerned by enrichment
        (a node is concerned by enrichment if one of his connected cell is enriched,
        ie if its evolution is modified by enrichment)
        """
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
        print(message)

    def compute_new_force(self, topologie, contrainte, classical_cell: np.array):
        """
        Calcul des forces agissant sur les noeuds
        :param topologie: topologie du calcul
        :param contrainte : tenseur des contriante de cauchy sigma xx
        :param classical_cell: masks of the classical cells
        :type topologie: Topology
        :type contrainte: numpy.array([nbr_of_node-1, 1], dtype=np.float64, order='C')
        """
        # Initialize node force to 0 because nodes force is now calculated with +=
        self._force = np.zeros_like(self._force)

        # For a node, force = stress on cell right - stress on cell left
        for cell_index in np.where(classical_cell)[0]:
            # Each cell gives a contribution to the force on node left and node right
            cell_connectivity = topologie.get_nodes_belonging_to_cell(cell_index)
            node_left = cell_connectivity[0]
            node_right = cell_connectivity[1]
            # if node is on left of the cell, cell is on right of the node
            if node_left != -1:
                self._force[node_left] += contrainte[cell_index] * self.section
            # if node is on right of the cell, cell is on left of the node
            if node_right != -1:
                self._force[node_right] -= contrainte[cell_index] * self.section

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
        self._upundemi[mask] = self._umundemi[mask] + \
                               multiplication_masse(matrice_masse, self.force[mask]) * delta_t

    def compute_complete_velocity_field(self):
        """
        Calcul du champ de vitesse vraie
        """
        self._v_field = np.copy(self._upundemi)

    def apply_correction_reference_bar(self, delta_t, inv_complete_mass_matrix,
                                       inv_wilkins_mass_matrix, mask):
        """
        Apply a correction on velocity field to compute velocity from exact(non lumped)
        mass matrix for elements in mask
        """
        self._upundemi[mask] += multiplication_masse(inv_complete_mass_matrix,
                                                     self.force[mask]) * delta_t
        self._upundemi[mask] -= multiplication_masse(inv_wilkins_mass_matrix,
                                                     self.force[mask]) * delta_t

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

