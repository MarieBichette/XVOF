#!/usr/bin/env python2.7
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

    def __init__(self, number_of_nodes, correction_3x3_on_cell_500=False):
        self.__number_of_nodes = number_of_nodes
        self.__mass_matrix = np.zeros([self.__number_of_nodes, 1], dtype=np.float64, order='C')
        self.__inv_mass_matrix = None
        self.__3x3_mass_matrix = np.zeros([3, 3])
        self.correction_3x3_on_cell_500 = correction_3x3_on_cell_500

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
        
    def compute_3x3_mass_matrix_for_cell_500(self, cell_mass_vector, mask=np.array([500])):
        """
        Compute the exact form of mass matrix (classical) , no lumping
        :param cell_mass_vector: vector of cells mass
        :param mask: id of cell to be considered
        """
        mass_498 = cell_mass_vector[mask-2]
        mass_499 = cell_mass_vector[mask-1]
        mass_500 = cell_mass_vector[mask]
        self.__3x3_mass_matrix[0, 0] = 3 * mass_498 + 2 * mass_499
        self.__3x3_mass_matrix[0, 1] = mass_499
        self.__3x3_mass_matrix[1, 0] = self.__3x3_mass_matrix[0, 1]
        self.__3x3_mass_matrix[1, 1] = 2 * mass_500 + 2*mass_499
        self.__3x3_mass_matrix[1, 2] = mass_500
        self.__3x3_mass_matrix[2, 1] = self.__3x3_mass_matrix[1, 2]
        self.__3x3_mass_matrix[2, 2] = 2 * mass_500
        self.__3x3_mass_matrix *= 1. / 6.
        self.__inv_3x3_mass_matrix = inverseMasse(self.__3x3_mass_matrix)

    @property
    def mass_matrix(self):
        """
        Accessor on _mass_matrix
        """
        return self.__mass_matrix

    @property
    def mass_matrix_3x3(self):
        """
        Accessor on _mass_matrix
        """
        return self.__3x3_mass_matrix

    @property
    def inverse_mass_matrix(self):
        """
        Accessor on the inverse of the mass matrix

        :return: the inverse of the mass matrix
        """
        return self.__inv_mass_matrix

    @property
    def inverse_3x3_mass_matrix(self):
        """
        Accessor on the inverse of the exact mass matrix

        :return: the inverse of the exact mass matrix
        """
        return self.__inv_3x3_mass_matrix
