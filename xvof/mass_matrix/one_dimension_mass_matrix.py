# -*- coding: iso-8859-1 -*-
"""
Implementing the OneDimensionMassMatrix class
"""

import numpy as np
from xvof.mass_matrix.mass_matrix import compute_wilkins_mass_matrix
from xvof.mass_matrix.mass_matrix_utilities import inverseMasse

class OneDimensionMassMatrix(object):
    """
    A class for 1d mass matrix
    """

    def __init__(self, number_of_nodes):
        self.__number_of_nodes = number_of_nodes
        self.__mass_matrix = np.zeros([self.__number_of_nodes, 1], dtype=np.float64, order='C')
        self.__inv_mass_matrix = None

    def compute_mass_matrix(self, topology, cell_mass_vector, node_number_by_cell_vector):
        """
        Compute the mass matrix and its inverse according to Wilkins method

        :param topology: topology of the simulation
        :param cell_mass_vector: cells mass vector
        :param node_number_by_cell_vector: number of nodes per cell (vector)
        """
        self.__mass_matrix = compute_wilkins_mass_matrix(topology, cell_mass_vector,
                                                         node_number_by_cell_vector)
        self.__inv_mass_matrix = inverseMasse(self.__mass_matrix)

    @property
    def inverse_mass_matrix(self):
        """
        Accessor on the inverse of the mass matrix

        :return: the inverse of the mass matrix
        """
        return self.__inv_mass_matrix

    @property
    def mass_matrix(self):
        """
        Accessor on _mass_matrix
        """
        return self.__mass_matrix


