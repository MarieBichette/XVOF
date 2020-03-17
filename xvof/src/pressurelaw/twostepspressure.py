#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant une pression à deux paliers constants
"""
from xvof.src.pressurelaw.pressurelaw import PressureLaw


class TwoStepsPressure(PressureLaw):
    """
    Une pression à deux paliers constants
    """
    def __init__(self, first_value, second_value, critical_time):
        super(TwoStepsPressure, self).__init__()
        self.__first_value = first_value
        self.__second_value = second_value
        self.__critical_time = critical_time

    def evaluate(self, time):
        if time < self.__critical_time:
            return self.__first_value
        else:
            return self.__second_value
