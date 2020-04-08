#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe g�rant la topologie du maillage
"""
import numpy as np


class Topology:
    """
    Une classe g�rant la topologie et les connectivit�s d'un maillage
    """
    def __init__(self, nbr_of_nodes, nbr_of_cells, dim=1):
        #
        self._dim = dim
        self._nbr_of_nodes = nbr_of_nodes
        self._nbr_of_cells = nbr_of_cells
        # Array dont chaque item est un array des indices des noeuds appartenant � la maille
        self._nodes_belonging_to_cell = np.ones(shape=(nbr_of_cells, 2 ** self._dim),
                                                dtype=np.int64, order='C') * (-1)
        # Array dont chaque item est un array des indices des mailles en contact avec le noeud
        self._cells_in_contact_with_node = np.ones(shape=(nbr_of_nodes, 2 ** self._dim),
                                                   dtype=np.int64, order='C') * (-1)

    @property
    def dimension(self):
        """
        Return the dimension of the mesh
        """
        return self._dim

    @property
    def nodes_belonging_to_cell(self):
        """
        Returns the array for nodes <-> cells connectivity
        """
        return self._nodes_belonging_to_cell

    @property
    def cells_in_contact_with_node(self):
        """
        Returns the array for cells <-> nodes connectivity
        """
        return self._cells_in_contact_with_node

    def set_nodes_belonging_to_cell(self, ind_cell, ind_node_list):
        """
        Attribue la liste des indices des noeuds, 'nodes_list'
        appartenant � la maille d'indice 'ind_cell'
        :param ind_cell: indice de la maille � laquelle attribuer les indices de noeuds
        : ind_node_list: indices des noeuds appartenant � la cellule d'indice ind_cell
        :type ind_cell: int
        :type ind_node_list: list
        """
        self._nodes_belonging_to_cell[ind_cell] = np.array(ind_node_list)

    def add_cell_in_contact_with_node(self, ind_node, ind_cell):
        """
        Ajoute l'indice, 'ind_cell', de la maille � la liste des mailles en contact
        avec le noeud d'indice 'ind_node'
        """
        k = 0
        while self._cells_in_contact_with_node[ind_node, k] != -1:
            k += 1
        self._cells_in_contact_with_node[ind_node, k] = ind_cell

    def get_nodes_belonging_to_cell(self, ind_cell):
        """
        Renvoie un tableau des noeuds appartenant � la maille
        :param ind_cell: indice de la maille dont on veut connaitre les noeuds
        :return: un tableau des noeuds appartenant � la maille
        :type ind_cell: int
        :rtype: numpy.array
        """
        return self._nodes_belonging_to_cell[ind_cell]

    def get_cells_in_contact_with_node(self, ind_node):
        """
        Renvoie un tableau des mailles en contact avec le noeud

        :param ind_node: indice du noeud dont on veut connaitre les mailles connexes
        :return: un tableau des mailles en contact avec le noeud

        :type ind_node: int
        :rtype: numpy.array
        """
        return self._cells_in_contact_with_node[ind_node]
