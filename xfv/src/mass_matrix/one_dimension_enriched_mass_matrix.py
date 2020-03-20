#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Implementing the OneDimensionEnrichedMassMatrix class
"""
import numpy as np

from xfv.src.mass_matrix.mass_matrix_utilities import SymNDArray, inverseMasse


class OneDimensionEnrichedMassMatrix(object):

    def __init__(self, size_of_mass_matrix):
        self._matrix_size = size_of_mass_matrix
        self._enriched_mass_matrix = \
            SymNDArray((size_of_mass_matrix, size_of_mass_matrix), dtype=np.float64, order='C')
        self._enriched_mass_matrix[:, :] = 0.
        self._inv_enriched_mass_matrix = \
            SymNDArray((size_of_mass_matrix, size_of_mass_matrix), dtype=np.float64, order='C')
        self._analytical_inv_enriched_mass_matrix = \
            SymNDArray((size_of_mass_matrix, size_of_mass_matrix), dtype=np.float64, order='C')

    def assemble_enriched_mass_matrix(self, *sub_matrix_names):
        """
        Assemble and inverse the mass matrix after enrichment
        """
        for name in sub_matrix_names:
            self._enriched_mass_matrix += getattr(self, name)
        print "Numerical inverse of mass matrix"
        self._inv_enriched_mass_matrix = inverseMasse(self._enriched_mass_matrix)

    @property
    def enriched_mass_matrix(self):
        """
        Accessor on the complete enriched_mass_matrix
        :return: the enriched mass matrix for a single discontinuity
        """
        return self._enriched_mass_matrix

    @property
    def inverse_enriched_mass_matrix(self):
        """
        Accessor on the inverse of the mass matrix
        :return: the inverse of the mass matrix
        """
        return self._inv_enriched_mass_matrix

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
        return self._inv_enriched_mass_matrix[self._matrix_size - 2:self._matrix_size, self._matrix_size - 2:self._matrix_size]

    @property
    def inverse_enriched_mass_matrix_coupling_dof(self):
        """
        Accessor on the inverse of the mass matrix for coupling between classical and enriched degrees of freedom
        :return: the coupling part of the inverse of the mass matrix
        """
        return self._inv_enriched_mass_matrix[0:self._matrix_size - 2, self._matrix_size - 2:self._matrix_size]


