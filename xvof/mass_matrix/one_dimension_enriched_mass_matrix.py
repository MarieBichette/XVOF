#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
A class defining the mass matrix 1D (enrichment)
"""
import numpy as np

from xvof.discontinuity.discontinuity import discontinuity_list
from xvof.mass_matrix.mass_matrix_utilities import SymNDArray, lump_matrix
from xvof.mass_matrix.mass_matrix import MassMatrix


class OneDimensionEnrichedMassMatrix(object):
    """
    Class for enriched mass matrix
    UO    Cell0    U1     Cell1      U2    Cell2    U3
    *------------ * ------ // ------ * ------------ *
    La matrice est organisée:
        U0    U1    U2    U3    U1 *  U2 *
    U0
    U1
    U2
    U3
    U1 *
    U2 *
    """
    def __init__(self, lumped_matrix_classic_dof=False, lumped_matrix_coupling=False, lumped_matrix_enr_dof=False):
        self._enriched_mass_matrix = SymNDArray((6,6), dtype=np.float64, order='C')
        self._enriched_mass_matrix[:, :] = 0.
        self._lumped_enriched_mass_matrix = SymNDArray((6,6), dtype=np.float64, order='C')
        self._lumped_enriched_mass_matrix[:, :] = 0.
        self._matrix_classic_dof = SymNDArray((6,6), dtype=np.float64, order='C')
        self._matrix_classic_dof_lumped = lumped_matrix_classic_dof
        self._matrix_coupling = SymNDArray((6,6), dtype=np.float64, order='C')
        self._matrix_coupling_lumped = lumped_matrix_coupling
        self._matrix_enr_dof = SymNDArray((6,6), dtype=np.float64, order='C')
        self._matrix_enr_dof_lumped = lumped_matrix_enr_dof

    def compute_enriched_mass_matrix_classical_dof(self, mass_moins_un, mass_0, mass_1, mass_2, mass_3):
        """
        Calcul la matrice de masse non condensée après enrichissement pour les degrés de libertés classiques
        Calcule la partie U0  U1  U2  U3
                        U0
                        U1
                        U2
                        U3
        """
        self._matrix_classic_dof[:,:] = 0.
        self._matrix_classic_dof[0, 0] = 2. * mass_0 + 6. * mass_moins_un / 2.
        self._matrix_classic_dof[0, 1] = mass_0
        self._matrix_classic_dof[1, 1] = 2* mass_0 + 2* mass_1
        self._matrix_classic_dof[1, 2] = mass_1
        self._matrix_classic_dof[2, 2] = 2 * mass_1 + 2 * mass_2
        self._matrix_classic_dof[2, 3] = mass_2
        self._matrix_classic_dof[3, 3] = 2 * mass_2 + 6. * mass_3 / 2.
        self._matrix_classic_dof *= 1. / 6.
        if self._matrix_classic_dof_lumped:
            self._matrix_classic_dof = lump_matrix(self._matrix_classic_dof)

    def compute_enriched_mass_matrix_enriched_dof(self, mass_0, mass_1, mass_2):
        """
        Calcul la matrice de masse non condensée après enrichissement pour les dégrés de liberté enrichis
        Calcule la partie U1* U2*
                        U1*
                        U2*
        """
        self._matrix_enr_dof[:, :] = 0.
        self._matrix_enr_dof[4, 4] = 2. * mass_0 + 2. * mass_1
        self._matrix_enr_dof[4, 5] = mass_1
        self._matrix_enr_dof[5, 5] = 2. * mass_1 + 2. * mass_2
        self._matrix_enr_dof *= 1. / 6.
        if self._matrix_enr_dof_lumped:
            self._matrix_enr_dof = lump_matrix(self._matrix_enr_dof)

    def compute_enriched_mass_matrix_couplage(self, mass_0, mass_1, mass_2, alpha):
        """
        Calcul la matrice de masse non condensée après enrichissement pour les couplages
        Calcule la partie U1* U2* et U0  U1  U2  U3
                       U0          U1*
                       U1          U2*
                       U2
                       U3
        """
        self._matrix_coupling[:, :] = 0.
        self._matrix_coupling[0, 4] = -mass_0
        self._matrix_coupling[1, 4] = -2 * mass_0 + mass_1 * (0.5 - 4  * alpha)
        self._matrix_coupling[1, 5] = (1 - 2 * alpha) * mass_1
        self._matrix_coupling[2, 4] = (1 - 2 * alpha) * mass_1
        self._matrix_coupling[2, 5] = (7. / 2. - 4. * alpha) * mass_1 + 2. * mass_2
        self._matrix_coupling[3, 5] = mass_2
        self._matrix_coupling *= 1. / 6.
        if self._matrix_coupling_lumped:
            self._matrix_coupling = lump_matrix(self._matrix_coupling)

    def compute_enriched_mass_matrix(self, topology, cells_mass):
        """
        Assemble la matrice de masse non condensée après enrichissement
        """
        # Boucle sur les discontinuités
        for disc in [d for d in discontinuity_list if not d.mass_matrix_updated]:
            print "Entre dans la boucle enriched mass pour la  discontinuite {:d}".format(disc.label)
            alpha = disc.position_in_ruptured_element
            # Suppose les éléments voisins triés par position croissante
            connectivity = topology.cells_in_contact_with_node[:]
            mask_in_nodes = disc.mask_in_nodes
            mask_out_nodes = disc.mask_out_nodes
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
            self.compute_enriched_mass_matrix_classical_dof(mass_moins_un, mass_0, mass_1, mass_2, mass_3)
            self.compute_enriched_mass_matrix_enriched_dof(mass_0, mass_1, mass_2)
            self.compute_enriched_mass_matrix_couplage(mass_0, mass_1, mass_2, alpha)

    def assemble_enriched_mass_matrix(self, *sub_matrix_names):
        """
        Assemblage de la matrice masse complète enrichie
        """
        for name  in sub_matrix_names:
            self._enriched_mass_matrix += getattr(self, name)

    @property
    def complete_mass_matrix(self):
        """
        Accessor on the complete enriched_mass_matrix
        """
        return self._enriched_mass_matrix

    def print_enriched_mass_matrix(self):
        m = self.complete_mass_matrix
        print "Enriched mass matrix :"
        ligne0 = "{:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}".format(m[0,0], m[0,1], m[0,2], m[0,3], m[0,4], m[0,5])
        ligne1 = "{:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}".format(m[1,0], m[1,1], m[1,2], m[1,3], m[1,4], m[1,5])
        ligne2 = "{:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}".format(m[2,0], m[2,1], m[2,2], m[2,3], m[2,4], m[2,5])
        ligne3 = "{:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}".format(m[3,0], m[3,1], m[3,2], m[3,3], m[3,4], m[3,5])
        ligne4 = "{:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}".format(m[4,0], m[4,1], m[4,2], m[4,3], m[4,4], m[4,5])
        ligne5 = "{:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}   {:+8.7g}".format(m[5,0], m[5,1], m[5,2], m[5,3], m[5,4], m[5,5])
        print ligne0
        print ligne1
        print ligne2
        print ligne3
        print ligne4
        print ligne5

