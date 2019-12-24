#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Implementation of MaximalStressCriterion class
"""
from xvof.rupturecriterion.rupturecriterion import RuptureCriterion


class MaximalStressCriterion(RuptureCriterion):
    """
    A rupture criterion based on minimal pressure
    """
    def __init__(self, sigma_max):
        super(MaximalStressCriterion, self).__init__()
        self.__max_stress = sigma_max

    def checkCriterion(self, cells, *args, **kwargs):
        """
        Return the mask of the cells where pressure is below the minimum pressure

        :param cells: cells on which to check the criterion
        :return: the mask of the cells where pressure is below the minimum pressure
        """
        return cells.stress[:, 0] > self.__max_stress


