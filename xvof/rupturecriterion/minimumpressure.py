# -*- coding: iso-8859-1 -*-
"""
Implementation of MinimumPressureCriterion class
"""
from xvof.rupturecriterion.rupturecriterion import RuptureCriterion


class MinimumPressureCriterion(RuptureCriterion):
    """
    A rupture criterion based on minimal pressure
    """
    def __init__(self, pmin):
        super(MinimumPressureCriterion, self).__init__()
        self.__minimum_pressure = pmin

    def checkCriterion(self, cells, *args, **kwargs):
        """
        Return the mask of the cells where pressure is below the minimum pressure

        :param cells: cells on which to check the criterion
        :return: the mask of the cells where pressure is below the minimum pressure
        """
        return cells.pressure.new_value < self.__minimum_pressure
