#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant un critère de rupture basé sur la pression minimum
"""
from xvof.rupturecriterion.rupturecriterion import RuptureCriterion


class MinimumPressureCriterion(RuptureCriterion):
    """
    Un critère de rupture basé sur la pression minimale
    """
    def __init__(self, pmin):
        RuptureCriterion.__init__(self)
        self.__minimum_pressure = pmin

    def checkCriterion(self, cell, *args, **kwargs):
        if cell.pression_t_plus_dt < self.__minimum_pressure:
            return True
        else:
            return False
