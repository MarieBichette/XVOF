#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe g�rant la topologie 1D du maillage
"""
import numpy as np

from xfv.src.mesh.topology import Topology


class Topology1D(Topology):
    """
    Sp�cialisation 1D de topology

    >>> my_topo = Topology1D(11, 10)
    >>> my_topo.get_cells_in_contact_with_node(0)
    array([ -1, 0])
    >>> my_topo.get_cells_in_contact_with_node(5)
    array([4, 5])
    >>> my_topo.get_nodes_belonging_to_cell(9)
    array([ 9, 10])
    >>> my_topo.get_cells_in_contact_with_node(np.array([1, 3]))
    array([[0, 1],
           [2, 3]])
    >>> my_topo.cells_in_contact_with_node[:]
    array([[ -1, 0],
           [ 0,  1],
           [ 1,  2],
           [ 2,  3],
           [ 3,  4],
           [ 4,  5],
           [ 5,  6],
           [ 6,  7],
           [ 7,  8],
           [ 8,  9],
           [ 9, -1]])
    """
    def __init__(self, nbr_of_nodes, nbr_of_cells):
        """
        Build a Class Topology1D for 1D meshes
        :param nbr_of_nodes: Number of nodes in the mesh
        :param nbr_of_cells: Number of cells in the mesh
        """
        Topology.__init__(self, nbr_of_nodes, nbr_of_cells)
        self._generate_mesh(nbr_of_cells)

    def add_cell_in_contact_with_node(self, ind_node, ind_cell):
        """
        Ajoute l'indice, 'ind_cell', de la maille � la liste des mailles en contact
        avec le noeud d'indice 'ind_node'
        """
        Topology.add_cell_in_contact_with_node(self, ind_node, ind_cell)
        for ind in range(self._nbr_of_nodes):
            conn = self._cells_in_contact_with_node[ind_node]
            if conn.size > 2:
                raise RuntimeError(f"The node {ind} is connected to more "
                                   f"than two cells {conn}.\n"
                                   "It is not correct in 1D!")

    def _generate_mesh(self, nbr_of_cells):
        """
        Generation du maillage (connectivity cells <-> nodes)
        """
        for ind_cell in range(nbr_of_cells):
            ind_node_left = ind_cell
            ind_node_right = ind_cell + 1
            self.add_cell_in_contact_with_node(ind_node_left, ind_cell)
            self.add_cell_in_contact_with_node(ind_node_right, ind_cell)
            self.set_nodes_belonging_to_cell(ind_cell, [ind_node_left, ind_node_right])

        # Avec cette construction, le premier noeud est connect� aux mailles [0, -1]
        # Or on voudrait que �a soit l'inverse pour le calcul des forces
        # => Swap :
        assert self._cells_in_contact_with_node[0, 1] == -1
        self._cells_in_contact_with_node[0, 1] = self._cells_in_contact_with_node[0, 0]
        self._cells_in_contact_with_node[0, 0] = -1
