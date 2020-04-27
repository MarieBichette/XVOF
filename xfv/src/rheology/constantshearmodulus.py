# -*- coding: utf-8 -*-
"""
Implementation of the ConstantShearModulus class
"""

from xfv.src.rheology.shearmodulus import ShearModulus


class ConstantShearModulus(ShearModulus):  # pylint: disable=too-few-public-methods
    """
    Class for constant shear modulus
    """

    def __init__(self, init_value):
        """
        Init of the class
        :param init_value: Value of the shear modulus
        """
        super(ConstantShearModulus, self).__init__(init_value)

    @classmethod
    def compute(cls):
        """
        Compute the shear modulus => returns constant value of shear modulus
        """
        return cls().shear_modulus
