# -*- coding: utf-8 -*-
"""
Implementing the EnrichedMassMatrixLump class
"""
import numpy as np
from abc import abstractmethod

from xfv.src.mass_matrix.enriched_mass_matrix import EnrichedMassMatrix
from xfv.src.mass_matrix.mass_matrix_utilities import inverse_masse


class EnrichedMassMatrixLump(EnrichedMassMatrix):
    """
    A class for the lumped enriched mass matrix
    """
    def __init__(self):
        """
        Build the class
        """
        matrix_size = 4
        super(EnrichedMassMatrixLump, self).__init__(matrix_size)
        self._enriched_mass_matrix = np.zeros([matrix_size], dtype=np.float64, order='C')
        self._inv_enriched_mass_matrix = np.zeros([matrix_size], dtype=np.float64, order='C')
        self._enriched_mass_matrix_left_part = np.zeros([matrix_size], dtype=np.float64, order='C')
        self._enriched_mass_matrix_right_part = np.zeros([matrix_size], dtype=np.float64, order='C')

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
    def inverse_enriched_mass_matrix_classic_dof(self) -> np.array:
        """
        Accessor on the inverse of the mass matrix for classical degrees of freedom
        :return: the extraction of the inverse of the mass matrix for classical dof
        """
        return self._inv_enriched_mass_matrix[0:self._matrix_size - 2]

    @property
    def inverse_enriched_mass_matrix_enriched_dof(self) -> np.array:
        """
        Accessor on the inverse of the mass matrix for enriched degrees of freedom
        :return: extraction of the inverse of the mass matrix for enriched dof
        """
        return self._inv_enriched_mass_matrix[self._matrix_size - 2:self._matrix_size]

    def get_mass_matrix_left(self) -> np.array:
        """
        Accessor on the part of mass matrix concerning the left part of the cracked cell
        """
        return self._enriched_mass_matrix_left_part

    def get_mass_matrix_right(self) -> np.array:
        """
        Accessor on the part of mass matrix concerning the right part of the cracked cell
        """
        return self._enriched_mass_matrix_right_part

    def rearrange_dof_in_inv_mass_matrix(self):
        """
        Rearrange dof to easily compute the node velocity with classical and enriched dof
        separately
        """
        res = np.zeros([self._matrix_size], dtype=np.float64, order='C')
        res[0] = self._inv_enriched_mass_matrix[0]
        res[1] = self._inv_enriched_mass_matrix[2]
        res[2] = self._inv_enriched_mass_matrix[1]
        res[3] = self._inv_enriched_mass_matrix[3]
        self._inv_enriched_mass_matrix = np.copy(res)

    def print_enriched_mass_matrix(self):
        """
        Print the mass matrix * (with aligned members)
        """
        m = self._enriched_mass_matrix  # pylint: disable=invalid-name
        print("Enriched mass matrix :")
        lign_0 = "{:+8.7g}   0.   0.   0.".format(m[0])
        lign_1 = "0.   {:+8.7g}   0.   0.".format(m[1])
        lign_2 = "0.   0.   {:+8.7g}   0.".format(m[2])
        lign_3 = "0.   0.   0.   {:+8.7g}".format(m[3])
        print(lign_0)
        print(lign_1)
        print(lign_2)
        print(lign_3)

    @abstractmethod
    def compute_enriched_mass_matrix_left_part(self, mass_0: float, mass_1: float, epsilon: float):
        """
        Compute the Hansbo mass matrix for the left part
        DDL are organized : 0 : N1g and 1 : N2g
        :param mass_0: mass of the element right on the left of the cracked cell
        :param mass_1 : mass of the cracked cell
        :param epsilon: relative position of the disc inside the cracked cell
        """
        pass

    @abstractmethod
    def compute_enriched_mass_matrix_right_part(self, mass_1: float, mass_2: float, epsilon: float):
        """
        Compute the Hansbo mass matrix for the right part
        DDL are organized : 2 : N2d and 3: N1d
        :param mass_1 : mass of the cracked cell
        :param mass_2: mass of the element right on the right of the cracked cell
        :param epsilon: relative position of the disc inside the cracked cell
        """
        pass
