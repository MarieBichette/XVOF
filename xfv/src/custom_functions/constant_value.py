# -*- coding: utf-8 -*-
"""
Implements the ConstantValue class
"""
from xfv.src.custom_functions.custom_function import CustomFunction


class ConstantValue(CustomFunction):
    """
    This class defines a function that returns a constant value

    .. image:: constant_value.png
        :scale: 75 %
        :align: center
    """
    def __init__(self, value):
        self.__value = value

    def evaluate(self, time, *args, **kwargs):
        """
        Returns the value of the function evaluated at time

        :param time: the required time
        :return: the value
        """
        return self.__value
