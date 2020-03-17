#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Implementation of MinimumPressureCriterion class
"""
from xvof.src.plasticitycriterion.plasticitycriterion import PlasticityCriterion
from xvof.src.utilities.stress_invariants_calculation import compute_J2


class VonMisesCriterion(PlasticityCriterion):
    """
    A plasticity criterion based on Von Mises model
    """
    def __init__(self):
        super(VonMisesCriterion, self).__init__()

    def checkCriterion(self, cells, *args, **kwargs):
        """
        Return the mask of the cells where pressure is below the minimum pressure

        :param cells: cells on which to check the criterion
        :return: the mask of the cells where plastic criterion is activated
        """
        J2 = compute_J2(cells.deviatoric_stress_new)
        return J2 > cells.yield_stress.new_value

    def checkCriterion_on_right_part_cells(self, disc, *args, **kwargs):
        """
        Return the mask of the cells where pressure is below the minimum pressure
        :param disc: current discontinuity
        :return: the mask of the cells where pressure is below the minimum pressure
        """
        J2 = compute_J2(disc.additional_dof_deviatoric_stress_new)
        return J2 > disc.additional_dof_yield_stress.new_value
