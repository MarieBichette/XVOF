# -*- coding: utf-8 -*-
"""
Implements the TwoSteps class
"""
from xfv.src.custom_functions.custom_function import CustomFunction


class TwoSteps(CustomFunction):
    """
    This class defines a 2 constant steps function

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
