#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe gérant la topologie 1D du maillage
"""
from xvof.mesh.topology import Topology


class Topology1D(Topology):
    '''
    Spécialisation 1D de topology
    '''
    def __init__(self, nbr_of_nodes, nbr_of_cells):
        Topology.__init__(self, nbr_of_nodes, nbr_of_cells)
        self._generateMesh(nbr_of_cells)

    def _generateMesh(self, nbr_of_cells):
        '''
        Generation du maillage (connectivité mailles <--> noeuds)
        '''
        for ind_cell in xrange(nbr_of_cells):
            ind_node_left = ind_cell
            ind_node_right = ind_cell + 1
            self.addCellInContactWithNode(ind_node_left, ind_cell)
            self.addCellInContactWithNode(ind_node_right, ind_cell)
            self.setNodesBelongingToCell(ind_cell, [ind_node_left, ind_node_right])

    def isRightBoundary(self, node):
        '''
        Renvoie vrai si le noeud 'node' est le noeud frontiere droite
        '''
        if node.index == len(self._nodes) - 1:
            return True
        else:
            return False

    def isLeftBoundary(self, node):
        '''
        Renvoie vrai le noeud 'node' est le noeud frontiere gauche
        '''
        if node.index == 0:
            return True
        else:
            return False
