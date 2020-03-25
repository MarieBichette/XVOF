#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Implementation d'une classe de module de cisaillement constant
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