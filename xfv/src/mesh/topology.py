# -*- coding: utf-8 -*-
"""
Module to manage mesh topology and objects connectivity
"""
import numpy as np


class Topology:
    """
    A class to manage mesh topology and objects connectivity
    """
    def __init__(self, nbr_of_nodes, nbr_of_cells, dim=1):
        #
        self._dim = dim
        self._nbr_of_nodes = nbr_of_nodes
        self._nbr_of_cells = nbr_of_cells
        # Array dont chaque item est un array des indices des noeuds appartenant ï¿½ la maille
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
        Register a node list as belonging to the cell 'ind_cell'

        :param ind_cell: cell index
        :param ind_node_list: list of the node index connected to ind_cell
        :type ind_cell: int
        :type ind_node_list: list
        """
        self._nodes_belonging_to_cell[ind_cell] = np.array(ind_node_list)

    def add_cell_in_contact_with_node(self, ind_node: int, ind_cell: int):
        """
        Register a cell ind_cell as belonging to the node ind_node

        :param ind_cell: cell index connected to ind_node
        :param ind_node: node index connected to ind_cell
        """
        conn = self._cells_in_contact_with_node[ind_node]
        first_emplace = np.argwhere(conn == -1)[0]
        conn[first_emplace[0]] = ind_cell

    def get_nodes_belonging_to_cell(self, ind_cell: int):
        """
        Returns the node indexes connected to the cell ind_cell
        :param ind_cell: cell_index
        """
        return self._nodes_belonging_to_cell[ind_cell]

    def get_cells_in_contact_with_node(self, ind_node: int):
        """
        Returns the cell indexes connected to the node ind_node

        :param ind_node: node index
        """
        return self._cells_in_contact_with_node[ind_node]
