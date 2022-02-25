#!/usr/bin/env python3.7
# -*- coding: iso-8859-1 -*-
"""
Implementing a non local criterion for failure + the cell is cracked if it verifies the failure criterion itself
"""
import numpy as np
from xfv.src.rupturecriterion.rupturecriterion import RuptureCriterion
from xfv.src.rupturecriterion.nonlocalstress import compute_neighbour


class NonLocalStressCriterionWithMax(RuptureCriterion):  # pylint: disable=too-few-public-methods
    """
    An class for non local failure criterion
    """
    def __init__(self, value, radius, average_strategy):
        super().__init__()
        self.critical_value = value
        self.radius = radius
        self.average_strategy= average_strategy

    def check_criterion(self, cells, *args, **kwargs):
        """
        Check of the rupture criterion on the cells in arguments
        :param cells: cells on which to check the criterion
        """
        mean_stress = np.zeros(cells.number_of_cells)
        nbr_div = np.zeros(cells.number_of_cells)

        for i in range(cells.number_of_cells):
            cells_in_radius, enr_cells_in_radius = compute_neighbour(cells, i, self.radius)

            mean_stress[i] += np.sum(cells.stress_xx[cells_in_radius])  \
                            + np.sum(cells.enr_stress_xx[enr_cells_in_radius])
            nbr_div[i] = len(np.where(cells_in_radius)[0]) + len(np.where(enr_cells_in_radius)[0])
        mean_stress = mean_stress / nbr_div
        return (mean_stress >= self.critical_value) * (cells.stress[:, 0] >= self.critical_value)
