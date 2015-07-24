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
    def __init__(self, list_of_nodes, list_of_cells):
        Topology.__init__(self, list_of_nodes, list_of_cells)
        self._generateMesh()

    def _generateMesh(self):
        '''
        Generation du maillage (connectivité mailles <--> noeuds)
        '''
        # On trie les noeuds selon les x croissants
        self._nodes = sorted(self._nodes, key=lambda m: m.coordt)
        # On affecte l'indice global des noeuds
        for ind, node in enumerate(self._nodes):
            node.index = ind
        #
        for ind, cell in enumerate(self._cells):
            cell.index = ind
            node_left = self._nodes[ind]
            node_right = self._nodes[ind + 1]
            self._addCellInContactWithNode(node_left, cell)
            self._addCellInContactWithNode(node_right, cell)
            self._setNodesBelongingToCell(cell, [node_left, node_right])

    def _isRightBoundary(self, node):
        '''
        Renvoie vrai si le noeud 'node' est le noeud frontiere droite
        '''
        if node.index == len(self._nodes) - 1:
            return True
        else:
            return False

    def _isLeftBoundary(self, node):
        '''
        Renvoie vrai le noeud 'node' est le noeud frontiere gauche
        '''
        if node.index == 0:
            return True
        else:
            return False
