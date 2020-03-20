#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Implementation d'une classe de limite d'élasticité (interface)
"""


class YieldStress(object):
    """
    Interface de limite d'élasticité
    """
    def __init__(self, initial_value, **kwargs):
        self.yield_stress = initial_value

    @classmethod
    def compute(cls):
        """
        Compute the new value of shear modulus
        :return : float
        """
        pass