# -*- coding: utf-8 -*-
"""
Implementation of the ConstantYieldStress class
"""

from xfv.src.rheology.yieldstress import YieldStress
import numpy as np


class ConstantYieldStress(YieldStress):  # pylint: disable=too-few-public-methods
    """
    A class for constant yield stress calculation
    """

    def __init__(self, init_value):
        """
        Initialization of the constant yield stress class
        :param init_value: initial yield stress
        """
        super(ConstantYieldStress, self).__init__(init_value)

    def compute(self, density: np.array) -> np.array:
        """
        Compute the value of the yield stress
        :return: float
        """
        return np.ones_like(density) * self.yield_stress
