#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe gérant la topologie du maillage
"""
from xvof.element.element1d import Element1d
from xvof.node.node import Node


class Topology(object):
    '''
    Une classe gérant la topologie et les connectivités d'un maillage
    '''
    def __init__(self, list_of_nodes):
        self._nodes = list_of_nodes[:]  # Liste des noeuds
        self._cells = []  # Liste des mailles
        # Liste dont chaque item est la liste des indices des noeuds appartenant à la maille
        self._nodes_belonging_to_cell = []
        # Liste dont chaque item est la liste des indices des mailles en contact avec le noeud
        self._cells_in_contact_with_node = []

    @property
    def cells(self):
        return self._cells

    @property
    def nodes(self):
        return self._nodes

    def _setNodesBelongingToCell(self, cell, nodes_list):
        '''
        Attribue la liste des indices des noeuds, 'nodes_list'
        appartenant à la maille 'cell'
        '''
        self._nodes_belonging_to_cell[cell.index] = [ nod.index for nod in nodes_list ]

    def _addNodeBelongingToCell(self, cell, node):
        '''
        Ajoute l'indice du noeud 'node' à la liste des noeuds
        appartenant à la maille 'cell'
        '''
        self._nodes_belonging_to_cell[cell.index].append(node.index)

    def _setCellsInContactWithNode(self, node, cells_list):
        '''
        Attribue la liste des indices des mailles, 'cells_list'
        en contact avec le noeud 'node'
        '''
        self._cells_in_contact_with_node[node.index] = [ cell.index for cell in cells_list ]

    def _addCellInContactWithNode(self, node, cell):
        '''
        Ajoute l'indice de la maille 'cell' à la liste des mailles en contact
        avec le noeud 'node'
        '''
        self._cells_in_contact_with_node[node.index].append(cell.index)

    def _getNodesBelongingToCell(self, cell):
        '''
        Renvoie la liste des noeuds appartenant à la maille
        '''
        res = [self._nodes[ind] for ind in self._nodes_belonging_to_cell[cell.index]]
        return res

    def _getCellsInContactWithNode(self, node):
        '''
        Renvoie la liste des mailles en contact avec le noeud
        '''
        res = [self._cells[ind] for ind in self._cells_in_contact_with_node[node.index]]
        return res


class Topology1D(Topology):
    '''
    Spécialisation 1D de topology
    '''
    def __init__(self, list_of_nodes, properties):
        Topology.__init__(self, list_of_nodes)
        self._cells_in_contact_with_node = [[] for _ in self._nodes]
        self._nodes_belonging_to_cell = [[] for _ in xrange(len(self._nodes) - 1)]
        self._generateMesh(properties)

    def _generateMesh(self, proprietes):
        # On trie les noeuds selon les x croissants
        self._nodes = sorted(self._nodes, key = lambda m : m.coordt)
        # On affecte l'indice global des noeuds
        for ind, node in enumerate(self._nodes):
            node.index = ind
        #
        stop = False
        ind_node = 0
        ind_cell = 0
        while not stop:
            node_left = self._nodes[ind_node]
            node_right = self._nodes[ind_node + 1]
            elem = Element1d(proprietes)
            elem.index = ind_cell
            self._cells.append(elem)
            self._addCellInContactWithNode(node_left, elem)
            self._addCellInContactWithNode(node_right, elem)
            self._setNodesBelongingToCell(elem, [node_left, node_right])
            if ind_node == len(self._nodes) - 2:
                stop = True
            ind_node += 1
            ind_cell += 1

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

