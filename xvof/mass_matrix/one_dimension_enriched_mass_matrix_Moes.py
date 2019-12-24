#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Implementing the OneDimensionEnrichedMassMatrix class for Moes enrichment
"""
import numpy as np

from xvof.mass_matrix.mass_matrix_utilities import SymNDArray
from xvof.mass_matrix.one_dimension_enriched_mass_matrix import OneDimensionEnrichedMassMatrix


class OneDimensionMoesEnrichedMassMatrix(OneDimensionEnrichedMassMatrix):
    """
    Class for 1D enriched mass matrix
    UO    Cell0    U1     Cell1      U2    Cell2    U3
    *------------ * ------ // ------ * ------------ *
    The matrix is organized as :
        U0    U1    U2    U3    U1 *  U2 *
    U0
    U1
    U2
    U3
    U1 *
    U2 *
    """
    def __init__(self, lump=None, correction_1_2=False):
        matrix_size = 6
        super(OneDimensionMoesEnrichedMassMatrix, self).__init__(matrix_size)
        # Mass matrix parts definition :
        self._matrix_classic_dof = SymNDArray((self._matrix_size, self._matrix_size), dtype=np.float64, order='C')
        self._matrix_coupling = SymNDArray((self._matrix_size, self._matrix_size), dtype=np.float64, order='C')
        self._matrix_enr_dof = SymNDArray((self._matrix_size, self._matrix_size), dtype=np.float64, order='C')
        self._analytical_inverse_of_mass_matrix = \
            SymNDArray((self._matrix_size, self._matrix_size), dtype=np.float64, order='C')

        self.lump = lump

        self.correction_1_2 = correction_1_2
        if self.correction_1_2:
            print " <!> APPLY CORRECTION 1/2 ON DIAGONAL TO RESPECT MASS CONSERVATION "

    def compute_enriched_mass_matrix_classical_dof(self, mass_moins_un, mass_0, mass_1, mass_2, mass_3):
        """
        Compute the mass matrix after enrichment for classical degrees of freedom
        Compute the part   U0  U1  U2  U3
                        U0
                        U1
                        U2
                        U3
        """
        self._matrix_classic_dof[:, :] = 0.
        if self.lump is not None:
            self._matrix_classic_dof[0, 0] = mass_0 / 2. + mass_moins_un / 2.
            self._matrix_classic_dof[1, 1] = mass_1 / 2. + mass_0 / 2.
            self._matrix_classic_dof[2, 2] = mass_2 / 2. + mass_1 / 2.
            self._matrix_classic_dof[3, 3] = mass_3 / 2. + mass_2 / 2.
            if self.correction_1_2:
                self._matrix_classic_dof[1, 1] = mass_1 / 4. + mass_0 / 2.
                self._matrix_classic_dof[2, 2] = mass_2 / 2. + mass_1 / 4.
        else:
            self._matrix_classic_dof[0, 0] = 2. * mass_0 + 3. * mass_moins_un
            self._matrix_classic_dof[0, 1] = mass_0
            self._matrix_classic_dof[1, 1] = 2. * mass_0 + 2. * mass_1
            self._matrix_classic_dof[1, 2] = mass_1
            self._matrix_classic_dof[2, 2] = 2. * mass_1 + 2. * mass_2
            self._matrix_classic_dof[2, 3] = mass_2
            self._matrix_classic_dof[3, 3] = 2. * mass_2 + 3. * mass_3
            self._matrix_classic_dof *= 1. / 6.

    def compute_enriched_mass_matrix_enriched_dof(self, mass_0, mass_1, mass_2):
        """
        Compute the mass matrix after enrichment for enriched degrees of freedom
        Compute the part    U1* U2*
                        U1*
                        U2*
        """
        self._matrix_enr_dof[:, :] = 0.
        if self.lump is not None:
            print "--> Pas d'assemblage sur les ddl enrichis"
            self._matrix_enr_dof[4, 4] = mass_1 / 2.
            self._matrix_enr_dof[5, 5] = mass_1 / 2.
            if self.correction_1_2:
                self._matrix_enr_dof[4, 4] = mass_1 / 4.
                self._matrix_enr_dof[5, 5] = mass_1 / 4.
        else:
            self._matrix_enr_dof[4, 4] = 2. * mass_0 + 2. * mass_1
            self._matrix_enr_dof[4, 5] = mass_1
            self._matrix_enr_dof[5, 5] = 2. * mass_1 + 2. * mass_2
            self._matrix_enr_dof *= 1. / 6.

    def compute_enriched_mass_matrix_couplage(self, mass_0, mass_1, mass_2, epsilon):
        """
        Compute the mass matrix after enrichment for coupled degrees of freedom
        Compute the part U1* U2* et U0  U1  U2  U3
                       U0          U1*
                       U1          U2*
                       U2
                       U3
        """
        alpha = 2 * epsilon-1
        self._matrix_coupling[:, :] = 0.
        if self.lump == "menouillard":
            self._matrix_coupling[1, 4] = mass_1 / 2. * (- 2 * epsilon + 1)
            self._matrix_coupling[2, 5] = mass_1 / 2. * (- 2 * epsilon + 1)
        # if self.lump == "diag_cst":
        #     # pas de termes de couplage
        #     pass
        if self.lump == "somme":
            self._matrix_coupling[1, 4] = mass_1 * (- 1 - 2 * alpha + alpha ** 2) / 4.
            self._matrix_coupling[2, 5] = mass_1 * (1 - 2 * alpha - alpha ** 2) / 4.

        else:
            self._matrix_coupling[0, 4] = -mass_0
            self._matrix_coupling[1, 4] = -2 * mass_0 + mass_1 * (2. - 12. * epsilon +
                                                                  12. * epsilon ** 2 - 4. * epsilon ** 3)
            self._matrix_coupling[1, 5] = (1. - 6. * epsilon ** 2 + 4. * epsilon ** 3) * mass_1
            self._matrix_coupling[2, 4] = (1. - 6. * epsilon ** 2 + 4. * epsilon ** 3) * mass_1
            self._matrix_coupling[2, 5] = (2. - 4. * epsilon ** 3) * mass_1 + 2. * mass_2
            self._matrix_coupling[3, 5] = mass_2
            self._matrix_coupling *= 1. / 6.

    def compute_enriched_mass_matrix(self, discontinuity, topology, cells_mass):
        """
        Compute the mass matrix after enrichment one discontinuity
        """
        # for disc in [d for d in Discontinuity.discontinuity_list() if not d.mass_matrix_updated]:
        print "Entre dans la boucle enriched mass pour la  discontinuite {:d}".format(discontinuity.label)
        print "Compute mass matrix with Moes method"
        epsilon = discontinuity.position_in_ruptured_element
        # Suppose les éléments voisins triés par position croissante
        connectivity = topology.cells_in_contact_with_node[:]
        mask_in_nodes = discontinuity.mask_in_nodes
        mask_out_nodes = discontinuity.mask_out_nodes
        cells_on_right = connectivity[mask_out_nodes][0]
        cells_on_left = connectivity[mask_in_nodes][0]
        cell_0 = cells_on_left[0]
        cell_1 = cells_on_left[1]
        if cell_1 != cells_on_right[0]:
            raise ValueError("Problème d'indice pour l'élément enrichi")
        cell_2 = cells_on_right[1]
        mass_moins_un = cells_mass[cell_0 - 1]
        mass_0 = cells_mass[cell_0]
        mass_1 = cells_mass[cell_1]
        mass_2 = cells_mass[cell_2]
        mass_3 = cells_mass[cell_2 + 1]
        self.compute_enriched_mass_matrix_classical_dof(mass_moins_un, mass_0, mass_1, mass_2, mass_3, epsilon)
        self.compute_enriched_mass_matrix_enriched_dof(mass_0, mass_1, mass_2, epsilon)
        self.compute_enriched_mass_matrix_couplage(mass_0, mass_1, mass_2, epsilon)

    def print_enriched_mass_matrix(self):
        """
        Print the mass matrix * (with aligned members)
        :return:
        """
        m = self.enriched_mass_matrix
        print "Enriched mass matrix :"
        ligne0 = "{:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}".format(m[0, 0], m[0, 1], m[0, 2],
                                                                                          m[0, 3], m[0, 4], m[0, 5])
        ligne1 = "{:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}".format(m[1, 0], m[1, 1], m[1, 2],
                                                                                          m[1, 3], m[1, 4], m[1, 5])
        ligne2 = "{:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}".format(m[2, 0], m[2, 1], m[2, 2],
                                                                                          m[2, 3], m[2, 4], m[2, 5])
        ligne3 = "{:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}".format(m[3, 0], m[3, 1], m[3, 2],
                                                                                          m[3, 3], m[3, 4], m[3, 5])
        ligne4 = "{:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}".format(m[4, 0], m[4, 1], m[4, 2],
                                                                                          m[4, 3], m[4, 4], m[4, 5])
        ligne5 = "{:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}".format(m[5, 0], m[5, 1], m[5, 2],
                                                                                          m[5, 3], m[5, 4], m[5, 5])
        print ligne0
        print ligne1
        print ligne2
        print ligne3
        print ligne4
        print ligne5
