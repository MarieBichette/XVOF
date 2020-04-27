# -*- coding: utf-8 -*-
"""
Implementation of the ConstantYieldStress class
"""

from xfv.src.rheology.yieldstress import YieldStress


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

    @classmethod
    def compute(cls):
        """
        Compute the value of the yield stress
        :return: float
        """
        return cls().yield_stress
