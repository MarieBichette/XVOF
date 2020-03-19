# -*- coding: utf-8 -*-
"""
Implements the ConstantVelocity class
"""
from xvof.src.boundary_condition.boundary_condition import BoundaryCondition


class ConstantVelocity(BoundaryCondition):
    """
    That class defines a constant velocity boundary condition
    """
    def __init__(self, value):
        self.__value = value

    def evaluate(self, time, *args, **kwargs):
        return self.__value
