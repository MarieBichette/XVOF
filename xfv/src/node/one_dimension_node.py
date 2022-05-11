# -*- coding: utf-8 -*-
"""
Module defining the OneDimensionNode class
"""
import numpy as np
from xfv.src.mass_matrix.mass_matrix_utilities import multiplication_masse
from xfv.src.data.data_container import DataContainer
from xfv.src.node import Node


class OneDimensionNode(Node):
    """
    A class to manage all the nodes in a 1d mesh
    """
    def __init__(self, nbr_of_nodes, poz_init, vit_init,
                 section=1.):

        super().__init__(nbr_of_nodes, position_initiale=poz_init,
                         dim=1, vitesse_initiale=vit_init)
        self._section = section
        self._nbr_of_nodes = nbr_of_nodes
        # By definition, these nodes are not enriched
        self._classical = np.empty([self.number_of_nodes, ], dtype=bool, order='C')
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
        By definition, a node is concerned by enrichment if one of his connected cell is enriched,
        ie if its evolution is modified by enrichment
        :return: boolean mask indicating which nodes are concerned by enrichment
        """
        return self._enrichment_not_concerned

    @property
    def velocity_field(self):
        """
        True field of the node velocities
        """
        return self._upundemi

    @property
    def section(self):
        """
        Surface associated with a node
        """
        return self._section

    def infos(self, index):
        """
        Affichage des informations concernant le noeud d'indice index

        :param index: node index to print
        :type index: int
        """
        Node.infos(self, index)
        message = "==> section = {:5.4g}".format(self.section)
        print(message)

    def compute_new_force(self, topologie, contrainte, classical_cell_mask: np.array):
        """
        Calcul des forces agissant sur les noeuds

        :param topologie: topologie du calcul
        :param contrainte: tenseur des contriante de cauchy sigma xx
        :param classical_cell_mask: masks of the classical cells
        :type topologie: Topology
        :type contrainte: numpy.array([nbr_of_node-1, 1], dtype=np.float64, order='C')
        """
        # Initialize node force to 0 because nodes force is now calculated with +=
        self._force = np.zeros_like(self._force)

        # Supposes that neighbors cells are sorted by increasing position
        node_in_contact_with_classical_cell = topologie.nodes_belonging_to_cell[classical_cell_mask]
        nodes_left = node_in_contact_with_classical_cell[:, 0]
        nodes_right = node_in_contact_with_classical_cell[:, 1]

        tmp_force = contrainte[classical_cell_mask] * self.section

        # For a node, force = stress on cell right - stress on cell left
        # if node is on left of the cell, cell is on right of the node
        self._force[:-1, 0][nodes_left] += tmp_force
        # if node is on right of the cell, cell is on left of the node
        self._force[0:, 0][nodes_right] -= tmp_force

    def compute_new_velocity(self, delta_t, mask, matrice_masse):
        """
        Computes the node velocities at time t+dt/2

        :param delta_t: pas de temps
        :param mask: boolean array to select non enriched nodes
        :param matrice_masse: inverse de la matrice de masse
        :type delta_t: float
        :type mask: boolean array

        """
        # ddl classique de noeud classique (sauf 0 1 2 3 quand enrichissement)
        # = noeuds classiques non concernï¿½s par l'enrichissement
        self._upundemi[mask] = self._umundemi[mask] + \
                               multiplication_masse(matrice_masse, self.force[mask]) * delta_t

    def compute_complete_velocity_field(self):
        """
        Compute the true field of nodal velocity
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
        Apply pressure on a given node

        :param ind_node: node index to apply the pressure at
        :param pressure: pressure to be applied
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
