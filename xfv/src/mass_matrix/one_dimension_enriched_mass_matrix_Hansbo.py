#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Implementing the OneDimensionEnrichedMassMatrix class
"""
import numpy as np

from xfv.src.mass_matrix.mass_matrix_utilities import SymNDArray
from xfv.src.mass_matrix.one_dimension_enriched_mass_matrix import OneDimensionEnrichedMassMatrix


class OneDimensionHansboEnrichedMassMatrix(OneDimensionEnrichedMassMatrix):
    def __init__(self, lump=None):
        matrix_size = 4
        super(OneDimensionHansboEnrichedMassMatrix, self).__init__(matrix_size)
        self._enriched_mass_matrix_left_part = \
            SymNDArray((self._matrix_size, self._matrix_size), dtype=np.float64, order='C')
        self._enriched_mass_matrix_left_part[:, :] = 0
        self._enriched_mass_matrix_right_part = \
            SymNDArray((self._matrix_size, self._matrix_size), dtype=np.float64, order='C')
        self._enriched_mass_matrix_right_part[:, :] = 0
        self.lump = lump

    def get_mass_matrix_left(self):
        return self._enriched_mass_matrix_left_part

    def get_mass_matrix_right(self):
        return self._enriched_mass_matrix_right_part

    def compute_enriched_mass_matrix_left_part(self, mass_0, mass_1, epsilon):
        """
        Compute the Hansbo mass matrix for the left part
        DDL are organized : 0 : N1g and 1 : N2g
        """
        # if self.lump == "diag_cst":
        #     self._enriched_mass_matrix_left_part[0, 0] = mass_1 / 2. + mass_0 / 2.
        #     self._enriched_mass_matrix_left_part[1, 1] = mass_1 / 2.
        if self.lump == "menouillard":
            self._enriched_mass_matrix_left_part[0, 0] = epsilon * mass_1 / 2. + mass_0 / 2.
            self._enriched_mass_matrix_left_part[1, 1] = epsilon * mass_1 / 2.
        elif self.lump == "somme":
            self._enriched_mass_matrix_left_part[0, 0] = mass_0 / 2. + epsilon * mass_1 / 2. * (2 - epsilon)
            self._enriched_mass_matrix_left_part[1, 1] = mass_1 / 2. * epsilon * epsilon
        else:
            self._enriched_mass_matrix_left_part[0, 0] = \
                epsilon * mass_1 - epsilon ** 2 * mass_1 + epsilon ** 3 / 3. * mass_1 + mass_0 / 2.
            self._enriched_mass_matrix_left_part[0, 1] = epsilon ** 2 * mass_1 / 2. - epsilon ** 3 / 3. * mass_1
            self._enriched_mass_matrix_left_part[1, 1] = epsilon ** 3 / 3. * mass_1

    def compute_enriched_mass_matrix_right_part(self, mass_1, mass_2,  epsilon):
        """
        Compute the Hansbo mass matrix for the right part
        DDL are organized : 2 : N2d and 3: N1d
        """
        # if self.lump == "diag_cst":
        #     self._enriched_mass_matrix_right_part[2, 2] = mass_1 / 2. + mass_2 / 2.
        #     self._enriched_mass_matrix_right_part[3, 3] = mass_1 / 2.
        if self.lump == "menouillard":
            self._enriched_mass_matrix_right_part[2, 2] = (1 - epsilon) * mass_1 / 2. + mass_2 / 2.
            self._enriched_mass_matrix_right_part[3, 3] = (1 - epsilon) * mass_1 / 2.
        elif self.lump == "somme":
            self._enriched_mass_matrix_right_part[2, 2] = mass_2 / 2. + (1 - epsilon) * mass_1 / 2. * (1 + epsilon)
            self._enriched_mass_matrix_right_part[3, 3] = mass_1 / 2. * (1 - epsilon) * (1 - epsilon)
        else:
            self._enriched_mass_matrix_right_part[2, 2] = 1. / 3. * mass_1 - epsilon ** 3 / 3. * mass_1 + mass_2 / 2.
            self._enriched_mass_matrix_right_part[2, 3] = 1. / 6. * mass_1 - \
                                                          epsilon ** 2 * mass_1 / 2. + epsilon ** 3 / 3. * mass_1
            self._enriched_mass_matrix_right_part[3, 3] = 1. / 3. * mass_1 - \
                                                          epsilon * mass_1 + epsilon ** 2 * mass_1 \
                                                          - epsilon ** 3 / 3. * mass_1

    def compute_enriched_mass_matrix(self, discontinuity, topology, cells_mass):
        """
        Compute the enriched mass matrix for Hansbo shape functions (associated with 1 discontinuity)
        :param discontinuity : discontinuity to be considered
        :param topology: topology = connectivity
        :param cells_mass: array of cells mass
        :return:
        """
        print "Entre dans la boucle enriched mass pour la  discontinuite {:d}".format(discontinuity.label)
        print "Compute mass matrix with Hansbo method"
        epsilon = discontinuity.position_in_ruptured_element
        # Suppose les éléments voisins triés par position croissante
        connectivity = topology.cells_in_contact_with_node[:]
        mask_in_nodes = discontinuity.mask_in_nodes
        mask_out_nodes = discontinuity.mask_out_nodes
        cells_on_right = connectivity[mask_out_nodes][0]
        cells_on_left = connectivity[mask_in_nodes][0]
        cell_0 = cells_on_left[0]
        cell_1 = cells_on_left[1]
        cell_2 = cells_on_right[1]
        if cell_1 != cells_on_right[0]:
            raise ValueError("Problème d'indice pour l'élément enrichi")
        mass_0 = cells_mass[cell_0]
        mass_1 = cells_mass[cell_1]
        mass_2 = cells_mass[cell_2]

        self.compute_enriched_mass_matrix_left_part(mass_0, mass_1, epsilon)
        self.compute_enriched_mass_matrix_right_part(mass_1, mass_2, epsilon)

    def rearrange_dof_in_inv_mass_matrix(self):

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
        m = self.enriched_mass_matrix
        print "Enriched mass matrix :"
        ligne0 = "{:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}".format(m[0, 0], m[0, 1], m[0, 2], m[0, 3])
        ligne1 = "{:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}".format(m[1, 0], m[1, 1], m[1, 2], m[1, 3])
        ligne2 = "{:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}".format(m[2, 0], m[2, 1], m[2, 2], m[2, 3])
        ligne3 = "{:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}".format(m[3, 0], m[3, 1], m[3, 2], m[3, 3])
        print ligne0
        print ligne1
        print ligne2
        print ligne3

