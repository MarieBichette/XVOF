#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Implementation d'une classe de module de cisaillement constant
"""

from xvof.src.rheology.shearmodulus import ShearModulus


class ConstantShearModulus(ShearModulus):

    def __init__(self, init_value):
        super(ConstantShearModulus, self).__init__(init_value)

    @classmethod
    def compute(cls):
        return cls().shear_modulus