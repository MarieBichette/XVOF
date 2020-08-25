# -*- coding: utf-8 -*-
"""
Implements the TwoSteps class
"""
from xfv.src.custom_functions.custom_function import CustomFunction


class TwoSteps(CustomFunction):
    """
    This class defines a 2 constant steps function

    .. image:: two_steps.png
        :scale: 75 %
        :align: center
    """
    def __init__(self, first_value, second_value, critical_time):
        self.__first_value = first_value
        self.__second_value = second_value
        self.__critical_time = critical_time

    def evaluate(self, time, *args, **kwargs):
        """
        Returns the value of the function evaluated at time

        :param time: the required time
        :return: the value
        """
        if time <= self.__critical_time:
            return self.__first_value
        return self.__second_value
