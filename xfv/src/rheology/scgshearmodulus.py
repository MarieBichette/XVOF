# -*- coding: utf-8 -*-
"""
Implementation of the SCGShearModulus class
"""

import numpy as np
from xfv.src.rheology.shearmodulus import ShearModulus


class SCGShearModulus(ShearModulus):  # pylint: disable=too-few-public-methods
    """
    Class for scg shear modulus
    """

    def __init__(self, init_value, init_density, Gp_prime):
        """
        Init of the class

        :param init_value: Value of the shear modulus
        """
        super().__init__(init_value, init_density, Gp_prime)
        self.init_value = init_value
        self.init_density = init_density
        self.Gp_prime = Gp_prime
	

    def compute(self, density: np.array, pressure) -> np.array:
        """
        Compute the shear modulus => returns SCG value of shear modulus

        :param density: the current density
        :return: the computed shear modulus
        """
        is_pressure_below_zero = pressure < 0.
        #is_pressure_ok = np.logical(is_pressure_below_zero)
        pressure[is_pressure_below_zero] = 0.
        density_part = np.power(density/self.init_density,(1./3.))
        G = np.ones_like(density) * self.init_value *(1. + self.Gp_prime/self.init_value * pressure / density_part)
        #print('G_max=',max(G))
        #print('G_min=',min(G))
        return G
