# -*- coding: utf-8 -*-
"""
Implementing the EnrichedMassMatrixLumpMenouillard class for the enriched mass matrix with
Menouillard lump
"""
from xfv.src.mass_matrix.enriched_mass_matrix_lump import EnrichedMassMatrixLump


class EnrichedMassMatrixLumpMenouillard(EnrichedMassMatrixLump):
    """
    A class for the lumped enriched mass matrix for Menouillard lumping
    """
    def __init__(self):
        """
        Build the class EnrichedMassMatrixLumpMenouillard
        """
        super(EnrichedMassMatrixLumpMenouillard, self).__init__()

    def compute_enriched_mass_matrix_left_part(self, mass_0: float, mass_1: float, epsilon: float):
        """
        Compute the Hansbo mass matrix for the left part
        DDL are organized : 0 : N1g and 1 : N2g
        :param mass_0: mass of the element right on the left of the cracked cell
        :param mass_1 : mass of the cracked cell
        :param epsilon: relative position of the disc inside the cracked cell
        """
        self._enriched_mass_matrix_left_part[0] = epsilon * mass_1 / 2. + mass_0 / 2.
        self._enriched_mass_matrix_left_part[1] = epsilon * mass_1 / 2.

    def compute_enriched_mass_matrix_right_part(self, mass_1: float, mass_2: float, epsilon: float):
        """
        Compute the Hansbo mass matrix for the right part
        DDL are organized : 2 : N2d and 3: N1d
        :param mass_1 : mass of the cracked cell
        :param mass_2: mass of the element right on the right of the cracked cell
        :param epsilon: relative position of the disc inside the cracked cell
        """
        self._enriched_mass_matrix_right_part[2] = (1 - epsilon) * mass_1 / 2. + mass_2 / 2.
        self._enriched_mass_matrix_right_part[3] = (1 - epsilon) * mass_1 / 2.
