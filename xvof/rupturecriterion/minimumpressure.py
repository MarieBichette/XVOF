#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant un critère de rupture basé sur la pression minimum
"""
import numpy as np
from xvof.rupturecriterion.rupturecriterion import RuptureCriterion


class MinimumPressureCriterion(RuptureCriterion):
    """
    Un critère de rupture basé sur la pression minimale
    """
    def __init__(self, pmin):
        super(MinimumPressureCriterion, self).__init__()
        self.__minimum_pressure = pmin

    def checkCriterion(self, cells, *args, **kwargs):
        # return cells.pressure.new_value < self.__minimum_pressure
        mask_milieu = np.ndarray(cells.pressure.new_value.shape, dtype=np.bool,order='C')
        mask_milieu[:] = False
        mask_milieu[500] = True
        return mask_milieu
