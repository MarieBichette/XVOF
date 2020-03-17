#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant une pression en rampe
"""
from xvof.boundary_condition.pressurelaw import PressureLaw


class RampPressure(PressureLaw):
    """
    Une pression à deux paliers constants
    """
    def __init__(self, first_value, second_value, start_time, end_time):
        """
        Ramp pressure
        :param first_value: intial pressure value
        :param second_value: pressure value to reach
        :param start_time: time for starting the ramp
        :param end_time: time for reaching second_value pressure
        """
        super(RampPressure, self).__init__()
        self.__first_value = first_value
        self.__second_value = second_value
        self.__start_time = start_time
        self.__end_time = end_time

    def evaluate(self, time, *args, **kwargs):
        if time <= self.__start_time:
            return self.__first_value
        elif time > self.__end_time:
            return self.__second_value
        else:  # time correspond à un instant dans la pente
            slope = (self.__second_value - self.__first_value) / (self.__end_time - self.__start_time)
            return slope * (time - self.__start_time) + self.__first_value
