# -*- coding: utf-8 -*-
"""
Implementation of the ConstantShearModulus class
"""

import numpy as np
from xfv.src.rheology.shearmodulus import ShearModulus


class ConstantShearModulus(ShearModulus):  # pylint: disable=too-few-public-methods
    """
    Class for constant shear modulus
    """

    def __init__(self, init_value, init_density, Gp_prime):
        """
        Init of the class

        :param init_value: Value of the shear modulus
        """
        super().__init__(init_value, init_density, Gp_prime)
        self.init_value = init_value
        self.init_density = init_density

    def compute(self, density: np.array, pressure) -> np.array:
        """
        Compute the shear modulus => returns constant value of shear modulus

        :param density: the current density
        :return: the computed shear modulus
        """
        
        return np.ones_like(density) * self.init_value
