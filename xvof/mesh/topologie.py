#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe gérant la topologie du maillage
"""


class Topology1D(object):
    '''
    Une classe gérant la topologie et les connectivités d'un maillage
    '''
    def __init__(self, nbr_cells, nbr_nodes):
        self.__nodes_belonging_to_cell = [[] for _ in range(nbr_cells)]
        self.__cells_in_contact_with_node = [[] for _ in range(nbr_nodes)]

    def setNodesBelongingToCell(self, cell, nodes_list):
        '''
        Attribue la liste des noeuds, 'nodes_list'
        appartenant à la maille  'cell'
        '''
        self.__nodes_belonging_to_cell[cell.indice] = nodes_list

    def addNodeBelongingToCell(self, cell, node):
        '''
        Ajoute le noeud 'node' à la liste des noeuds 
        appartenant à la maille 'cell'
        '''
        self.__nodes_belonging_to_cell[cell.indice].append(node)

    def setCellsInContactWithNode(self, node, cells_list):
        '''
        Attribue la liste des mailles, 'cells_list'
        en contact avec le noeud 'node'
        '''
        self.__cells_in_contact_with_node[node.index] = cells_list

    def addCellInContactWithNode(self, node, cell):
        '''
        Ajoute la maille 'cell' à la liste des mailles en contact 
        avec le noeud 'node'
        '''
        self.__cells_in_contact_with_node[node.index].append(cell)

    def getNodesBelongingToCell(self, cell):
        '''
        Renvoie la liste des noeuds appartenant à la maille
        '''
        return self.__nodes_belonging_to_cell[cell.indice]

    def getCellsInContactWithNode(self, node):
        '''
        Renvoie la liste des mailles en contact avec le noeud
        '''
        return self.__cells_in_contact_with_node[node.index]
