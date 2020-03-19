# -*- coding: utf-8 -*-
"""
Implements the ConstantPressure class
"""
from xvof.src.boundary_condition.boundary_condition import BoundaryCondition


class ConstantPressure(BoundaryCondition):
    """
    This class defines a constant pressure boundary condition

         ^
    value|.................
         |
         |
         |
         |_________________>

    """
    def __init__(self, value):
        self.__value = value

    def evaluate(self, time, *args, **kwargs):
        return self.__value
