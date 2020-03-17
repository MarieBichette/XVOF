#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe gérant la topologie 1D du maillage
"""
import numpy as np

from xvof.mesh.topology import Topology


class Topology1D(Topology):
    """
    Spécialisation 1D de topology

    >>> my_topo = Topology1D(11, 10)
    >>> my_topo.getCellsInContactWithNode(0)
    array([ 0, -1])
    >>> my_topo.getCellsInContactWithNode(5)
    array([4, 5])
    >>> my_topo.getNodesBelongingToCell(9)
    array([ 9, 10])
    >>> my_topo.getCellsInContactWithNode(np.array([1, 3]))
    array([[0, 1],
           [2, 3]])
    >>> my_topo.cells_in_contact_with_node[:]
    array([[ 0, -1],
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
        Topology.__init__(self, nbr_of_nodes, nbr_of_cells)
        self._nbr_of_nodes = nbr_of_nodes
        self._nodes_belonging_to_cell = np.ndarray(shape=(nbr_of_cells, 2), dtype=np.int64, order='C')
        self._nodes_belonging_to_cell[:, :] = -1
        self._cells_in_contact_with_node = np.ndarray(shape=(nbr_of_nodes, 2), dtype=np.int64, order='C')
        self._cells_in_contact_with_node[:, :] = -1
        self._generateMesh(nbr_of_cells)

    def addCellInContactWithNode(self, ind_node, ind_cell):
        """
        Ajoute l'indice, 'ind_cell', de la maille à la liste des mailles en contact
        avec le noeud d'indice 'ind_node'
        """
        Topology.addCellInContactWithNode(self,ind_node, ind_cell)
        for ind in range(self._nbr_of_nodes):
            cells_in_contact_with_node_number = 0
            cells_in_contact_with_node_number += np.size(self.getCellsInContactWithNode(ind))
            if cells_in_contact_with_node_number > 2:
                raise SystemExit

    def _generateMesh(self, nbr_of_cells):
        """
        Generation du maillage (connectivité mailles <--> noeuds)
        """
        for ind_cell in xrange(nbr_of_cells):
            ind_node_left = ind_cell
            ind_node_right = ind_cell + 1
            self.addCellInContactWithNode(ind_node_left, ind_cell)
            self.addCellInContactWithNode(ind_node_right, ind_cell)
            self.setNodesBelongingToCell(ind_cell, [ind_node_left, ind_node_right])
