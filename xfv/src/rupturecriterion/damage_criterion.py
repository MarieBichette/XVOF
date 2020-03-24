#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Implementation of DamageCriterion class
"""
from xfv.src.rupturecriterion.rupturecriterion import RuptureCriterion


class DamageCriterion(RuptureCriterion):   # pylint: disable=too-few-public-methods
    """
    A rupture criterion based on damage value
    """
    def __init__(self, d_limite):
        super(DamageCriterion, self).__init__()
        self.__limit_damage = d_limite

    def check_criterion(self, cells, *args, **kwargs):
        """
        Return the mask of the cells where pressure is below the minimum pressure

        :param cells: cells on which to check the criterion
        :return: the mask of the cells where pressure is below the minimum pressure
        """
        return cells.damage_variable > self.__limit_damage
