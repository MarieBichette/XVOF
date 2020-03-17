#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Implementation d'une classe de limite d'élasticité constant
"""

from xvof.rheology.yieldstress import YieldStress


class ConstantYieldStress(YieldStress):

    def __init__(self, init_value):
        super(ConstantYieldStress, self).__init__(init_value)

    @classmethod
    def compute(cls):
        return cls().yield_stress