#!/usr/bin/env python3.7
# -*- coding: iso-8859-1 -*-
"""
Implementing a non local criterion for failure
"""
import numpy as np
from xfv.src.mass_matrix.mass_matrix_utilities import SymNDArray
from xfv.src.rupturecriterion.rupturecriterion import RuptureCriterion


def compute_weight(cells, weight_strategy):
    """
    Return an array containing the weight (between 0 and 1)
    Index [i, j] : weight for cell j with respect to cell i
    """
    mask = cells.cell_in_target
    coord = cells.coordinates_x[mask]
    size = len(coord)  # nombre de mailles dans la cible
    weight_matrix = SymNDArray((size, size), dtype=np.float64)  # matrice symétrique par construction

    enr_coord = np.copy(coord)
    mask_enriched = cells.enriched[mask]
    enr_coord[mask_enriched] = cells.enr_coordinates_x[mask_enriched]
    enr_weight_matrix = SymNDArray((size, size), dtype=np.float64)  # matrice symétrique par construction

    # On remplit le triangle superieur de la matrice, et par construction, le triangle inférieur se remplit aussi
    # Le calcul enrichi = même calcul sauf que l'on prend les coordonnées enrichies au lieu des coordonnées classiques
    # (si la maille est classique, on prendra deux fois la partie classique en compte)
    for i in range(1, size):
        weight_matrix[:i, i] = weight_strategy.compute_weight(np.abs(coord[:i] - coord[i]))
        enr_weight_matrix[:i, i] = weight_strategy.compute_weight(np.abs(enr_coord[:i] - enr_coord[i]))
        weight_matrix[i, i] = 1.
        enr_weight_matrix[i, i] = 1.
    # Reste le poids de la maille 0 par rapport à elle même...
    weight_matrix[0, 0] = 1.
    enr_weight_matrix[0, 0] = 1.

    return weight_matrix, enr_weight_matrix


class NonLocalStressCriterion(RuptureCriterion):  # pylint: disable=too-few-public-methods
    """
    An class for non local failure criterion
    """
    def __init__(self, value, radius, average_strategy):
        super().__init__()
        self.critical_value = value
        self.weight_strategy = average_strategy
        self.weight_strategy.set_radius(radius)

    def check_criterion(self, cells, *args, **kwargs):
        """
        Check of the rupture criterion on the cells in arguments
        :param cells: cells on which to check the criterion
        """
        weight_matrix, enr_weight_matrix = compute_weight(cells, self.weight_strategy)

        stress = cells.stress_xx[cells.cell_in_target]
        enr_stress = cells.enr_stress_xx[cells.cell_in_target]

        mean_stress = np.dot(weight_matrix, stress)
        enr_mean_stress = np.dot(enr_weight_matrix, enr_stress)
        nbr_div = np.sum(weight_matrix, axis=0) + np.sum(enr_weight_matrix, axis=0)

        mean_stress = (mean_stress + enr_mean_stress) / nbr_div
        return mean_stress >= self.critical_value
