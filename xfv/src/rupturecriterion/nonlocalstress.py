#!/usr/bin/env python3.7
# -*- coding: iso-8859-1 -*-
"""
Implementing a non local criterion for failure
"""
import numpy as np
from xfv.src.rupturecriterion.rupturecriterion import RuptureCriterion


class NonLocalStressCriterion(RuptureCriterion):  # pylint: disable=too-few-public-methods
    """
    An class for non local failure criterion
    """
    def __init__(self, value):
        super().__init__()
        self.critical_value = value
        self.stencile = 1

    def check_criterion(self, cells, *args, **kwargs):
        """
        Check of the rupture criterion on the cells in arguments
        :param cells: cells on which to check the criterion
        """
        mean_stress = np.zeros(cells.number_of_cells)  # taille nb_mailles
        # note : pas de verification du critere de rupture sur les mailles au bord de la geometrie
        for i in range(1, cells.number_of_cells - 1):  # 1 à nb_mailles - 1
            # import ipdb ; ipdb.set_trace()
            if cells.enriched[i]:
                stress_g = cells.enr_stress_xx[i]
            else:
                stress_g = cells.stress_xx[i - 1]
            stress = cells.stress_xx[i]
            stress_d = cells.stress_xx[i+1]
            mean_stress[i] = (stress_g + stress + stress_d) / 3.
        return mean_stress > self.critical_value