#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe gérant la topologie du maillage
"""
import numpy as np

class Topology(object):
    '''
    Une classe gérant la topologie et les connectivités d'un maillage
    '''
    def __init__(self, nbr_of_nodes, nbr_of_cells, dim=1):
        #
        self._dim = dim
        # Array dont chaque item est un array des indices des noeuds appartenant à la maille
        self._nodes_belonging_to_cell = [None for _ in xrange(nbr_of_cells)]
        # Array dont chaque item est un array des indices des mailles en contact avec le noeud
        self._cells_in_contact_with_node = [None for _ in xrange(nbr_of_nodes)]

    @property
    def dimension(self):
        return self._dim

    def setNodesBelongingToCell(self, ind_cell, ind_node_list):
        '''
        Attribue la liste des indices des noeuds, 'nodes_list'
        appartenant à la maille d'indice 'ind_cell'

        :param ind_cell: indice de la maille à laquelle attribuer les indices de noeuds
        : ind_node_list: indices des noeuds appartenant à la cellule d'indice ind_cell

        :type ind_cell: int
        :type ind_node_list: list
        '''
        self._nodes_belonging_to_cell[ind_cell] = np.array(ind_node_list)

    def addNodeBelongingToCell(self, ind_cell, ind_node):
        '''
        Ajoute l'indice du noeud 'node' à la liste des indices des noeuds
        appartenant à la maille d'indice 'ind_cell'

        :param ind_cell: indice de la maille à laquelle attribuer l'indice du noeud
        :param ind_node: indice du noeud à attribuer à la maille
        '''
        if self._nodes_belonging_to_cell[ind_cell] is not None:
            self._nodes_belonging_to_cell[ind_cell] = np.append(self._nodes_belonging_to_cell[ind_cell], ind_node)
        else:
            self._nodes_belonging_to_cell[ind_cell] = np.array([ind_node])

    def setCellsInContactWithNode(self, ind_node, ind_cell_list):
        '''
        Attribue la liste des indices des mailles, 'ind_cell_list'
        au noeud d'indice 'ind_node'

        :param ind_node: indice du noeud auquel attribuer les mailles
        :param ind_cell_list: liste des indices des mailles à attribuer au noeud

        :type ind_node: int
        :type ind_cell_list: list
        '''
        self._cells_in_contact_with_node[ind_node] = np.array(ind_cell_list)

    def addCellInContactWithNode(self, ind_node, ind_cell):
        '''
        Ajoute l'indice, 'ind_cell', de la maille à la liste des mailles en contact
        avec le noeud d'indice 'ind_node'
        '''
        if self._cells_in_contact_with_node[ind_node] is not None:
            self._cells_in_contact_with_node[ind_node] = np.append(self._cells_in_contact_with_node[ind_node], ind_cell)
        else:
            self._cells_in_contact_with_node[ind_node] = np.array([ind_cell])

    def getNodesBelongingToCell(self, ind_cell):
        '''
        Renvoie un tableau des noeuds appartenant à la maille

        :param ind_cell: indice de la maille dont on veut connaitre les noeuds
        :return: un tableau des noeuds appartenant à la maille

        :type ind_cell: int
        :rtype: numpy.array
        '''
        return self._nodes_belonging_to_cell[ind_cell]

    def getCellsInContactWithNode(self, ind_node):
        '''
        Renvoie un tableau des mailles en contact avec le noeud

        :param ind_node: indice du noeud dont on veut connaitre les mailles connexes
        :return: un tableau des mailles en contact avec le noeud

        :type ind_node: int
        :rtype: numpy.array
        '''
        return self._cells_in_contact_with_node[ind_node]

