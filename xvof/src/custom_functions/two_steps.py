# -*- coding: utf-8 -*-
"""
Implements the TwoStepsPressure class
"""
from xvof.src.boundary_condition.boundary_condition import BoundaryCondition


class TwoStepsPressure(BoundaryCondition):
    """
    This class defines a 2 constant steps pressure boundary condition

                ^
    second value|    ...........
                |    |
                |    |
    first value |....|
                |_________________>
                 critical time
    """
    def __init__(self, first_value, second_value, critical_time):
        self.__first_value = first_value
        self.__second_value = second_value
        self.__critical_time = critical_time

    def evaluate(self, time, *args, **kwargs):
        if time <= self.__critical_time:
            return self.__first_value
        return self.__second_value
