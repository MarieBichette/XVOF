# -*- coding: utf-8 -*-
"""
Implements the ConstantPressure class
"""
from xvof.src.boundary_condition.pressure_law import PressureLaw


class ConstantPressure(PressureLaw):
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
        super(ConstantPressure, self).__init__()
        self.__value = value

    def evaluate(self, time, *args, **kwargs):
        return self.__value
