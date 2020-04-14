#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""
Implementing the OneDimensionEnrichedMassMatrix class with total lump
"""
from xfv.src.mass_matrix.enriched_mass_matrix_lump import EnrichedMassMatrixLump


class EnrichedMassMatrixLumpSum(EnrichedMassMatrixLump):
    """
    A class for the enriched mass matrix lumped with sum of lines of the consistent enriched
    mass matrix
    """
    def __init__(self):
        """
        Build the class EnrichedMassMatrixLumpSum
        """
        super(EnrichedMassMatrixLumpSum, self).__init__()

    def compute_enriched_mass_matrix_left_part(self, mass_0: float, mass_1: float, epsilon: float):
        """
        Compute the Hansbo mass matrix for the left part
        DDL are organized : 0 : N1g and 1 : N2g
        :param mass_0: mass of the element right on the left of the cracked cell
        :param mass_1 : mass of the cracked cell
        :param epsilon: relative position of the disc inside the cracked cell
        """
        self._enriched_mass_matrix_left_part[0] = mass_0 / 2. \
                                                  + epsilon * mass_1 / 2. * (2 - epsilon)
        self._enriched_mass_matrix_left_part[1] = mass_1 / 2. * epsilon * epsilon

    def compute_enriched_mass_matrix_right_part(self, mass_1: float, mass_2: float, epsilon: float):
        """
        Compute the Hansbo mass matrix for the right part
        DDL are organized : 2 : N2d and 3: N1d
        :param mass_1 : mass of the cracked cell
        :param mass_2: mass of the element right on the right of the cracked cell
        :param epsilon: relative position of the disc inside the cracked cell
        """
        self._enriched_mass_matrix_right_part[2] = mass_2 / 2. \
                                                   + (1 - epsilon) * mass_1 / 2. * (1 + epsilon)
        self._enriched_mass_matrix_right_part[3] = mass_1 / 2. * (1 - epsilon) * (1 - epsilon)
