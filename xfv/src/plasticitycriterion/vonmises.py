# -*- coding: utf-8 -*-
"""
Implementation of MinimumPressureCriterion class
"""
from xfv.src.plasticitycriterion.plasticitycriterion import PlasticityCriterion
from xfv.src.utilities.stress_invariants_calculation import compute_second_invariant


class VonMisesCriterion(PlasticityCriterion):
    """
    A plasticity criterion based on Von Mises model
    """
    def check_criterion(self, cells):
        """
        Return the mask of the cells where the VonMises plasticity criterion is verified

        :param cells: cells on which to check the criterion
        :return: the mask of the cells where the VonMises plasticity criterion is verified
        """
        return (compute_second_invariant(cells.deviatoric_stress_new) >
                cells.yield_stress.new_value)

    @staticmethod
    def check_criterion_on_right_part_cells(cells):
        """
        Return the mask of the cells where the VonMises plasticity criterion is verified

        :param cells: cells on which to check the criterion
        :return: the mask of the cells where the VonMises plasticity criterion is verified
        """
        return (compute_second_invariant(cells.additional_dof_deviatoric_stress_new) >
                cells.additional_dof_yield_stress.new_value)
