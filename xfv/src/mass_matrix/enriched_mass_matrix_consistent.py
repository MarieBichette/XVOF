# -*- coding: utf-8 -*-
"""
Implementing the EnrichedMassMatrixConsistent class
"""
import numpy as np

from xfv.src.mass_matrix.mass_matrix_utilities import SymNDArray, inverse_masse
from xfv.src.mass_matrix.enriched_mass_matrix import EnrichedMassMatrix


class EnrichedMassMatrixConsistent(EnrichedMassMatrix):
    """
    A class for the consistent enriched mass matrix
    """
    def __init__(self):
        """
        Build the class
        """
        matrix_size = 4
        super(EnrichedMassMatrix, self).__init__(matrix_size)
        # Matrix for the enriched cell
        self._enriched_mass_matrix = \
            SymNDArray((self._matrix_size, self._matrix_size), dtype=np.float64, order='C')
        self._enriched_mass_matrix[:, :] = 0.
        # And its inverse
        self._inv_enriched_mass_matrix = \
            SymNDArray((self._matrix_size, self._matrix_size), dtype=np.float64, order='C')
        self._inv_enriched_mass_matrix[:, :] = 0.
        # Submatrix for the left part of the enriched cell
        self._enriched_mass_matrix_left_part = \
            SymNDArray((self._matrix_size, self._matrix_size), dtype=np.float64, order='C')
        self._enriched_mass_matrix_left_part[:, :] = 0
        # Submatrix for the right part of the enriched cell
        self._enriched_mass_matrix_right_part = \
            SymNDArray((self._matrix_size, self._matrix_size), dtype=np.float64, order='C')
        self._enriched_mass_matrix_right_part[:, :] = 0

    def get_mass_matrix_left(self):
        """
        Accessor on the part of mass matrix concerning the left part of the cracked cell
        """
        return self._enriched_mass_matrix_left_part

    def get_mass_matrix_right(self):
        """
        Accessor on the part of mass matrix concerning the right part of the cracked cell
        """
        return self._enriched_mass_matrix_right_part

    def assemble_enriched_mass_matrix(self, *sub_matrix_names):
        """
        Assemble and inverse the mass matrix after enrichment
        """
        for name in sub_matrix_names:
            self._enriched_mass_matrix += getattr(self, name)
        self.print_enriched_mass_matrix()
        print("Numerical inverse of mass matrix")
        self._inv_enriched_mass_matrix = inverse_masse(self._enriched_mass_matrix)

    @property
    def inverse_enriched_mass_matrix_classic_dof(self):
        """
        Accessor on the inverse of the mass matrix for classical degrees of freedom
        :return: the extraction of the inverse of the mass matrix for classical dof
        """
        return self._inv_enriched_mass_matrix[0:self._matrix_size - 2, 0:self._matrix_size - 2]

    @property
    def inverse_enriched_mass_matrix_enriched_dof(self):
        """
        Accessor on the inverse of the mass matrix for enriched degrees of freedom
        :return: extraction of the inverse of the mass matrix for enriched dof
        """
        return self._inv_enriched_mass_matrix[self._matrix_size - 2:self._matrix_size,
                                              self._matrix_size - 2:self._matrix_size]

    def compute_enriched_mass_matrix_left_part(self, mass_0: float, mass_1: float, epsilon: float):
        """
        Compute the Hansbo mass matrix for the left part
        DDL are organized : 0 : N1g and 1 : N2g
        :param mass_0: mass of the element right on the left of the cracked cell
        :param mass_1 : mass of the cracked cell
        :param epsilon: relative position of the disc inside the cracked cell
        """
        self._enriched_mass_matrix_left_part[0, 0] = \
            epsilon * mass_1 - epsilon ** 2 * mass_1 + epsilon ** 3 / 3. * mass_1 + mass_0 / 2.
        self._enriched_mass_matrix_left_part[0, 1] = \
            epsilon ** 2 * mass_1 / 2. - epsilon ** 3 / 3. * mass_1
        self._enriched_mass_matrix_left_part[1, 1] = epsilon ** 3 / 3. * mass_1

    def compute_enriched_mass_matrix_right_part(self, mass_1: float, mass_2: float, epsilon: float):
        """
        Compute the Hansbo mass matrix for the right part
        DDL are organized : 2 : N2d and 3: N1d
        :param mass_1 : mass of the cracked cell
        :param mass_2: mass of the element right on the right of the cracked cell
        :param epsilon: relative position of the disc inside the cracked cell
        """
        self._enriched_mass_matrix_right_part[2, 2] = \
            1. / 3. * mass_1 - epsilon ** 3 / 3. * mass_1 + mass_2 / 2.
        self._enriched_mass_matrix_right_part[2, 3] = \
            1. / 6. * mass_1 - epsilon ** 2 * mass_1 / 2. + epsilon ** 3 / 3. * mass_1
        self._enriched_mass_matrix_right_part[3, 3] = \
            1. / 3. * mass_1 - epsilon * mass_1 + epsilon ** 2 * mass_1 - epsilon ** 3 / 3. * mass_1

    def rearrange_dof_in_inv_mass_matrix(self):
        """
        Rearrange dof to easily compute the node velocity with classical and enriched dof separately
        """
        res = SymNDArray((self._matrix_size, self._matrix_size), dtype=np.float64, order='C')
        res[:, :] = 0
        res[0, 0] = self._inv_enriched_mass_matrix[0, 0]
        res[1, 1] = self._inv_enriched_mass_matrix[2, 2]
        res[2, 2] = self._inv_enriched_mass_matrix[1, 1]
        res[3, 3] = self._inv_enriched_mass_matrix[3, 3]
        res[0, 2] = self._inv_enriched_mass_matrix[0, 1]
        res[1, 3] = self._inv_enriched_mass_matrix[2, 3]
        self._inv_enriched_mass_matrix = np.copy(res)

    def print_enriched_mass_matrix(self):
        """
        Print the mass matrix * (with aligned members)
        :return:
        """
        m = self._enriched_mass_matrix  # pylint: disable=invalid-name
        print("Enriched mass matrix :")
        lign_0 = "{:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}".format(
            m[0, 0], m[0, 1], m[0, 2], m[0, 3])
        lign_1 = "{:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}".format(
            m[1, 0], m[1, 1], m[1, 2], m[1, 3])
        lign_2 = "{:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}".format(
            m[2, 0], m[2, 1], m[2, 2], m[2, 3])
        lign_3 = "{:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}".format(
            m[3, 0], m[3, 1], m[3, 2], m[3, 3])
        print(lign_0)
        print(lign_1)
        print(lign_2)
        print(lign_3)
