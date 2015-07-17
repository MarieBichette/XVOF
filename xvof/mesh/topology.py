#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe gérant la topologie du maillage
"""


class Topology(object):
    '''
    Une classe gérant la topologie et les connectivités d'un maillage
    '''
    def __init__(self, list_of_nodes, list_of_cells):
        self._nodes = list_of_nodes[:]  # Liste des noeuds
        self._cells = list_of_cells[:]  # Liste des mailles
        # Liste dont chaque item est la liste des indices des noeuds appartenant à la maille
        self._nodes_belonging_to_cell = [[] for _ in self._cells]
        # Liste dont chaque item est la liste des indices des mailles en contact avec le noeud
        self._cells_in_contact_with_node = [[] for _ in self._nodes]

    @property
    def cells(self):
        '''
        Retourne une copie (superficielle) de la liste des mailles
        '''
        return self._cells[:]

    @property
    def nodes(self):
        '''
        Retourne une copie (superficielle) de la liste des noeuds
        '''
        return self._nodes[:]

    def _setNodesBelongingToCell(self, cell, nodes_list):
        '''
        Attribue la liste des indices des noeuds, 'nodes_list'
        appartenant à la maille 'cell'
        '''
        self._nodes_belonging_to_cell[cell.index] = [nod.index for nod in nodes_list]

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
        self._cells_in_contact_with_node[node.index] = [cell.index for cell in cells_list]

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

    def _changeNodeType(self, node):
        '''
        Change le type de noeud. Le noeud d'indice 'node.index' est remplacé par 'node'
        '''
        self._nodes[node.index] = node

    def _changeCellType(self, cell):
        '''
        Change le type de maille. La maille d'indice 'cell.index' est remplacée par 'cell'
        '''
        self._cells[cell.index] = cell
