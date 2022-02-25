#!/usr/bin/env python3.7
# -*- coding: iso-8859-1 -*-
"""
Implementing a non local criterion for failure
"""
import numpy as np
from xfv.src.rupturecriterion.rupturecriterion import RuptureCriterion


def compute_neighbour(cells, i, radius):
    """
    Return a mask of cells located in the neighbourhood of the target i
    """
    distance_cell_to_i, enr_distance_cell_to_i = compute_distance(cells, i)
    cells_in_radius = (distance_cell_to_i < radius).flatten()
    enr_cells_in_radius = (enr_distance_cell_to_i < radius).flatten()
    # Only enriched cells count in the computation of the "enriched" stress mean
    enr_cells_in_radius = enr_cells_in_radius * cells.enriched
    return cells_in_radius, enr_cells_in_radius


def compute_distance(cells, i):
    """
    Return the distance between cells and target i
    """
    distance_cell_to_i = np.abs(cells.coordinates_x - cells.coordinates_x[i])
    enr_distance_cell_to_i = np.abs(cells.coordinates_x - cells.enr_coordinates_x[i])
    return distance_cell_to_i, enr_distance_cell_to_i


class NonLocalStressCriterion(RuptureCriterion):  # pylint: disable=too-few-public-methods
    """
    An class for non local failure criterion
    """
    def __init__(self, value, radius, average_strategy):
        super().__init__()
        self.critical_value = value
        self.radius = radius
        self.average_strategy = average_strategy

    def check_criterion(self, cells, *args, **kwargs):
        """
        Check of the rupture criterion on the cells in arguments
        :param cells: cells on which to check the criterion
        """
        mean_stress = np.zeros(cells.number_of_cells)
        nbr_div = np.zeros(cells.number_of_cells)

        for i in range(cells.number_of_cells):
            distance_to_cell_i, enr_distance_to_cell_i = compute_distance(cells, i)
            weight_i = self.average_strategy.compute_weight(self.radius, distance_to_cell_i)
            enr_weight_i = self.average_strategy.compute_weight(self.radius, enr_distance_to_cell_i)
            # weight_i et #enr_weight_i are 0 for cells outside the neighbourhood
            mean_stress[i] = np.sum(cells.stress_xx * weight_i + cells.enr_stress_xx * enr_weight_i)
            nbr_div[i] = np.sum(weight_i + enr_weight_i)
        mean_stress = mean_stress / nbr_div
        return mean_stress >= self.critical_value
