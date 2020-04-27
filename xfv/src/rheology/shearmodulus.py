# -*- coding: utf-8 -*-
"""
Interface for the shear modulus calculation
"""


class ShearModulus:  # pylint: disable=too-few-public-methods
    """
    Abstract class for the shear modulus computation
    """
    def __init__(self, initial_value):
        self.shear_modulus = initial_value

    @classmethod
    def compute(cls):
        """
        Compute the new value of shear modulus
        """
