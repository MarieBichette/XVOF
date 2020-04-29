# -*- coding: utf-8 -*-
"""
Implementing the EnrichedMassMatrix class
"""
from abc import abstractmethod


class EnrichedMassMatrix(object):
    """
    A class to factorize code for the enriched mass matrix
    """

    def __init__(self, size_of_mass_matrix):
        """
        Build the class EnrichedMassMatrix
        :param size_of_mass_matrix: size of the matrix
        """
        self._matrix_size = size_of_mass_matrix

    def compute_enriched_mass_matrix(self, discontinuity, topology, cells_mass):
        """
        Compute the enriched mass matrix for Hansbo shape functions
        (associated with 1 discontinuity)
        :param discontinuity : discontinuity to be considered
        :param topology: topology = connectivity
        :param cells_mass: array of cells mass
        """
        print("Compute enriched mass matrix for discontinuity {:d}".format(discontinuity.label))
        epsilon = discontinuity.position_in_ruptured_element
        # Suppose les éléments voisins triés par position croissante
        connectivity = topology.cells_in_contact_with_node[:]
        cells_on_right = connectivity[discontinuity.mask_out_nodes][0]
        cells_on_left = connectivity[discontinuity.mask_in_nodes][0]
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

    @abstractmethod
    def get_mass_matrix_left(self):
        """
        Accessor on the part of mass matrix concerning the left part of the cracked cell
        """
        pass

    @abstractmethod
    def get_mass_matrix_right(self):
        """
        Accessor on the part of mass matrix concerning the right part of the cracked cell
        """
        pass

    @abstractmethod
    def assemble_enriched_mass_matrix(self, *sub_matrix_names):
        """
        Assemble and inverse the mass matrix after enrichment
        """
        pass

    @abstractmethod
    def inverse_enriched_mass_matrix_classic_dof(self):
        """
        Accessor on the inverse of the mass matrix for classical degrees of freedom
        :return: the extraction of the inverse of the mass matrix for classical dof
        """
        return

    @abstractmethod
    def inverse_enriched_mass_matrix_enriched_dof(self):
        """
        Accessor on the inverse of the mass matrix for enriched degrees of freedom
        :return: extraction of the inverse of the mass matrix for enriched dof
        """
        return

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

    @abstractmethod
    def rearrange_dof_in_inv_mass_matrix(self):
        """
        Rearrange dof to easily compute the node velocity with classical and enriched dof separately
        """
        pass

    @abstractmethod
    def print_enriched_mass_matrix(self):
        """
        Print the mass matrix * (with aligned members)
        :return:
        """
        pass
