# -*- coding: utf-8 -*-
"""
Implementation of the ConstantYieldStress class
"""

import numpy as np
from xfv.src.rheology.yieldstress import YieldStress


class SCGYieldStress(YieldStress):  # pylint: disable=too-few-public-methods
    """
    A class for constant yield stress calculation
    """

    def __init__(self, init_value, init_shear_modulus, Y_max, beta, m):
        """
        Initialization of the constant yield stress class

        :param init_value: initial yield stress
        """
        super().__init__(init_value, init_shear_modulus, Y_max, beta, m)
        self.init_value = init_value
        self.init_shear_modulus = init_shear_modulus
        self.Y_max = Y_max
        self.beta = beta
        self.m = m 

    def compute(self, density: np.array, strain_plastic_eq, G) -> np.array:
        """
        Compute the value of the yield stress

        :param density: the current density
        :return: the computed yield stress
        """
        #print('max epsilon_eq=',max(strain_plastic_eq))
        part_Y = np.ones_like(density) * self.init_value*(1.+self.beta*np.power(strain_plastic_eq,self.m))
        is_Ymax_below_partY = self.Y_max < part_Y
        #is_Y_ok = np.logical(is_Ymax_below_partY)
        part_Y[is_Ymax_below_partY] = self.Y_max
        Y = part_Y*G/self.init_shear_modulus
        #print('Y_max',max(Y))
        #print('Y_min',min(Y))
        return Y
