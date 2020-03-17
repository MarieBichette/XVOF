#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant une pression en rampe
"""
from xvof.src.boundary_condition.pressurelaw import PressureLaw
from xvof.src.boundary_condition.ramppressure import RampPressure


class CreneauRampPressure(PressureLaw):
    """
    Une pression à deux paliers constants
    """
    def __init__(self, first_value, second_value, third_value, start_time, end_first_slope_time,
                 begin_second_slope_time, end_time):
        """
        Ramp pressure
        :param first_value: initial pressure value
        :param second_value: pressure value to reach on top of creneau
        :param third_value : pressure at the end
        :param start_time: time for starting the ramp
        :param end_first_slope_time : time to reach the top of creneau
        :param begin_second_slope_time : time to start second ramp
        :param end_time: time for reaching third_value pressure
        """
        super(CreneauRampPressure, self).__init__()
        self.__first_value = first_value
        self.__second_value = second_value
        self.__third_value = third_value
        self.__start_time = start_time
        self.__end_first_slope_time = end_first_slope_time
        self.__begin_second_slope_time = begin_second_slope_time
        self.__end_time = end_time

        if self.__end_first_slope_time < self.__start_time or \
                        self.__begin_second_slope_time < self.__end_first_slope_time or \
                        self.__end_time < self.__begin_second_slope_time:
            raise ValueError ("""Cannot go into the past. You have to build the ramp with increasing times""")

        self.__ramp_1 = RampPressure(self.__first_value, self.__second_value,
                                     self.__start_time, self.__end_first_slope_time)
        self.__ramp_2 = RampPressure(self.__second_value, self.__third_value,
                                     self.__begin_second_slope_time, self.__end_time)

    def evaluate(self, time, *args, **kwargs):
        if time <= self.__begin_second_slope_time:
            return self.__ramp_1.evaluate(time)
        else:
            return self.__ramp_2.evaluate(time)
