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
        self._nodes_belonging_to_cell = None
        # Array dont chaque item est un array des indices des mailles en contact avec le noeud
        self._cells_in_contact_with_node = None

    @property
    def dimension(self):
        return self._dim

    @property
    def nodes_belonging_to_cell(self):
        return self._nodes_belonging_to_cell

    @property
    def cells_in_contact_with_node(self):
        return self._cells_in_contact_with_node

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

    def addCellInContactWithNode(self, ind_node, ind_cell):
        '''
        Ajoute l'indice, 'ind_cell', de la maille à la liste des mailles en contact
        avec le noeud d'indice 'ind_node'
        '''
        k = 0
        while self._cells_in_contact_with_node[ind_node, k] != -1:
            k += 1
        self._cells_in_contact_with_node[ind_node, k] = ind_cell

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

