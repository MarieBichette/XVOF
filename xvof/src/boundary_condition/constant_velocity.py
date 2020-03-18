# -*- coding: utf-8 -*-
"""
Implements the ConstantVelocity class
"""
from xvof.src.boundary_condition.velocity_law import VelocityLaw


class ConstantVelocity(VelocityLaw):
    """
    That class defines a constant velocity boundary condition
    """
    def __init__(self, value):
        super(ConstantVelocity, self).__init__()
        self.__value = value

    def evaluate(self, time, *args, **kwargs):
        return self.__value
