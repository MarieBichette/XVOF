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
    size = len(coord)
    weight_matrix = np.zeros([size, size], dtype=np.float64)
    enr_weight_matrix = np.zeros([size, size], dtype=np.float64)

    # Le calcul enrichi = m�me calcul sauf que l'on prend les coordonn�es enrichies au lieu des coordonn�es classiques
    # (si la maille est classique, on prendra deux fois la partie classique en compte)
    for i in range(0, size):
        weight_matrix[:, i] = weight_strategy.compute_weight(np.abs(coord - coord[i])).flatten()

    if cells.enriched.any():
        enr_coord = np.copy(cells.coordinates_x)
        mask_enriched = np.where(np.logical_and(mask, cells.enriched))[0]
        enr_coord[mask_enriched] = cells.enr_coordinates_x[mask_enriched]
        for i in range(0, size):
            enr_weight_matrix[:, i] = weight_strategy.compute_weight(np.abs(enr_coord[cells.cell_in_target] - enr_coord[i])).flatten()

    return weight_matrix, enr_weight_matrix


class NonLocalStressCriterion(RuptureCriterion):  # pylint: disable=too-few-public-methods
    """
    An class for non local failure criterion
    """
    def __init__(self, value, average_strategy):
        super().__init__()
        self.critical_value = value
        self.weight_strategy = average_strategy

    def check_criterion(self, cells, *args, **kwargs):
        """
        Check of the rupture criterion on the cells in arguments
        :param cells: cells on which to check the criterion
        """
        # Pour gagner un peu de temps de calcul, on ne calcule le crit�re que sur la cible (pas sur le projectile)
        weight_matrix, enr_weight_matrix = compute_weight(cells, self.weight_strategy)

        stress = cells.stress_xx[cells.cell_in_target]
        enr_stress = np.copy(stress)
        mask_enriched = np.where(cells.enriched[cells.cell_in_target])[0]  # mailles enrichies de target
        enr_stress[mask_enriched] = cells.enr_stress_xx[cells.cell_in_target][mask_enriched]

        mean_stress = np.dot(weight_matrix, stress)
        enr_mean_stress = np.dot(enr_weight_matrix, enr_stress)
        nbr_div = np.sum(weight_matrix, axis=0) + np.sum(enr_weight_matrix, axis=0)
        # axis = 0 ou 1 (la matrice est sym�trique)
        global_stress = (mean_stress.flatten() + enr_mean_stress.flatten()) / nbr_div

        # On reprend la taille nb_cell
        return_value = np.zeros(cells.number_of_cells, dtype=bool)
        return_value[cells.cell_in_target] = global_stress >= self.critical_value
        return return_value
