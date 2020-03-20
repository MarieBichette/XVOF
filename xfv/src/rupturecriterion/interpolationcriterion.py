#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Implementation of InterpolationPressureCriterion class
"""
from xfv.src.rupturecriterion.rupturecriterion import RuptureCriterion


class InterpolationPressureCriterion(RuptureCriterion):
    """
    A rupture criterion based on minimal pressure
    """
    def __init__(self, pmin, stencil):
        super(InterpolationPressureCriterion, self).__init__()
        self.__minimum_pressure = pmin

    def checkCriterion(self, cells, *args, **kwargs):
        """
        Return the mask of the cells where pressure is below the minimum pressure

        :param cells: cells on which to check the criterion
        :return: the mask of the cells where pressure is below the minimum pressure
        """
        return cells.pressure.new_value < self.__minimum_pressure


