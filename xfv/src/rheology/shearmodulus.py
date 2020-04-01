#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Implementation d'une classe de module de cisaillement (interface)
"""


class ShearModulus(object):  # pylint: disable=too-few-public-methods
    """
    Classe m�re de module de cisaillement
    """
    def __init__(self, initial_value):
        self.shear_modulus = initial_value

    @classmethod
    def compute(cls):
        """
        Compute the new value of shear modulus
        :return : float
        """
        pass
