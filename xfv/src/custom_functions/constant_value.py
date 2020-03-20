# -*- coding: utf-8 -*-
"""
Implements the ConstantValue class
"""
from xfv.src.custom_functions.custom_function import CustomFunction


class ConstantValue(CustomFunction):
    """
    This class defines a function that returns a constant value

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
