#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Implementation of MaximalStressCriterion class
"""
from xfv.src.rupturecriterion.rupturecriterion import RuptureCriterion


class MaximalStressCriterion(RuptureCriterion):   # pylint: disable=too-few-public-methods
    """
    A rupture criterion based on minimal pressure
    """
    def __init__(self, sigma_max):
        super().__init__()
        self.__max_stress = sigma_max

    def check_criterion(self, cells, *args, **kwargs):
        """
        Return the mask of the cells where pressure is below the minimum pressure

        :param cells: cells on which to check the criterion
        :return: the mask of the cells where pressure is below the minimum pressure
        """
        return cells.stress[:, 0] > self.__max_stress
