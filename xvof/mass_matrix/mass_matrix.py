#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
A class defining the mass matrix
"""
import numpy as np

from abc import abstractmethod

from xvof.mass_matrix.mass_matrix_utilities import inverseMasse

class MassMatrix(object):
    """
    Class for mass matrix (general)
    """
    def __init__(self,number_of_nodes):
        self.number_of_nodes = number_of_nodes
        # import ipdb ; ipdb.set_trace()
        self._mass_matrix = np.zeros([self.number_of_nodes, 1], dtype=np.float64, order='C')


    def compute_mass_matrix(self, topologie, vecteur_masse_elements, vecteur_nb_noeuds_par_element):
        """
        Calcule les masses associées à chaque noeud par moyenne arithmétique de la
        masse des éléments voisins (méthode Wilkins)
        self.mass_matix is type array 1D (= vecteur)

        :param topologie: topologie du calcul
        :param vecteur_masse_elements: vecteur des masses de chaque élément
        :param vecteur_nb_noeuds_par_element: vecteur des nombres de noeuds que possède chaque élément

        :type topologie: Topology
        :type vecteur_masse_elements: numpy.array([nbr_of_nodes, 1], dtype=np.float64, order='C')
        :type vecteur_nb_noeuds_par_element: numpy.array([nbr_of_nodes, 1], dtype=np.int64, order='C')
        """
        for ind_node in xrange(self.number_of_nodes):
            elements_voisins = topologie.getCellsInContactWithNode(ind_node)
            # Indice -1 => element inexistant
            elements_voisins = elements_voisins[elements_voisins != -1]
            self._mass_matrix[ind_node] = np.sum(vecteur_masse_elements[elements_voisins] /
                                           vecteur_nb_noeuds_par_element[elements_voisins])
        return self._mass_matrix


    @property
    def mass_matrix_value(self):
        """Accessor on _mass_matrix"""
        return self._mass_matrix





