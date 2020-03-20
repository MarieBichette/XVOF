# -*- coding: iso-8859-1 -*-
"""
Implementig the HalfRodComparisonCriterion class
"""
import numpy as np
from xfv.src.rupturecriterion.rupturecriterion import RuptureCriterion


class HalfRodComparisonCriterion(RuptureCriterion):
    """
    A not physical rupture criterion to validate XFEM by comparing
    results to other results obtained without XFEM on a half rod
    """
    def __init__(self, ruptured_cell_index=500):
        super(HalfRodComparisonCriterion, self).__init__()
        self.__ruptured_cell_index = int(ruptured_cell_index)

    def checkCriterion(self, cells):
        """
        Return the mask of the cells where only the ruptured_cell_index is set to True

        :param cells: cells on which to check criterion
        :return: mask of the cells where only the ruptured_cell_index is set to True
        """
        mask_milieu = np.ndarray(cells.pressure.new_value.shape, dtype=np.bool, order='C')
        mask_milieu[:] = False
        mask_milieu[self.__ruptured_cell_index] = True
        return mask_milieu
