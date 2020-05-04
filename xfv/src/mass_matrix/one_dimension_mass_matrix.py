# -*- coding: utf-8 -*-
"""
Implementing the OneDimensionMassMatrix class
"""

import numpy as np
from xfv.src.mass_matrix.mass_matrix import compute_wilkins_mass_matrix
from xfv.src.mass_matrix.mass_matrix_utilities import inverse_masse


class OneDimensionMassMatrix:
    """
    A class for 1d mass matrix
    """

    def __init__(self, number_of_nodes, consistent_matrix_on_last_cells=False):
        self.__number_of_nodes = number_of_nodes
        self.__mass_matrix = np.zeros([self.__number_of_nodes, 1], dtype=np.float64, order='C')
        self.__inv_mass_matrix = None
        self.__correction_mass_matrix = None
        self.__inv_correction_mass_matrix = None
        self.consistent_mass_matrix_on_last_cells = consistent_matrix_on_last_cells

    def compute_mass_matrix(self, topology, cell_mass_vector, node_number_by_cell_vector):
        """
        Compute the mass matrix and its inverse according to Wilkins method

        :param topology: topology of the simulation
        :param cell_mass_vector: cells mass vector
        :param node_number_by_cell_vector: number of nodes per cell (vector)
        """
        self.__mass_matrix = compute_wilkins_mass_matrix(topology, cell_mass_vector,
                                                         node_number_by_cell_vector)
        self.__inv_mass_matrix = inverse_masse(self.__mass_matrix)

    def compute_correction_mass_matrix_for_cell_500(self, cell_mass_vector, mask_node, topologie):
        """
        Compute the exact form of mass matrix (classical) , no lumping
        :param cell_mass_vector: vector of cells mass
        :param mask_node: id of cell to be considered
        :param topologie : topology
        """
        connect = topologie.cells_in_contact_with_node[mask_node]
        mask_cell = np.unique(connect)[1:]
        shape = len(mask_cell)
        self.__correction_mass_matrix = np.zeros([shape, shape])
        mass_499 = cell_mass_vector[mask_cell][-2]
        mass_500 = cell_mass_vector[mask_cell][-1]
        self.__correction_mass_matrix[0, 0] = 2 * mass_500 + 3 * mass_499
        self.__correction_mass_matrix[0, 1] = mass_500
        self.__correction_mass_matrix[1, 0] = self.__correction_mass_matrix[0, 1]
        self.__correction_mass_matrix[1, 1] = 2 * mass_500
        self.__correction_mass_matrix *= 1. / 6.
        self.__inv_correction_mass_matrix = inverse_masse(self.__correction_mass_matrix)

    @property
    def mass_matrix(self):
        """
        Accessor on _mass_matrix
        """
        return self.__mass_matrix

    @property
    def mass_matrix_correction(self):
        """
        Accessor on _mass_matrix
        """
        return self.__correction_mass_matrix

    @property
    def inverse_mass_matrix(self):
        """
        Accessor on the inverse of the mass matrix

        :return: the inverse of the mass matrix
        """
        return self.__inv_mass_matrix

    @property
    def inverse_correction_mass_matrix(self):
        """
        Accessor on the inverse of the exact mass matrix

        :return: the inverse of the exact mass matrix
        """
        return self.__inv_correction_mass_matrix
