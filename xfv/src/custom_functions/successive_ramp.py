# -*- coding: utf-8 -*-
"""
Implements the SuccessiveRamp class
"""
from xfv.src.custom_functions.custom_function import CustomFunction


class SuccessiveRamp(CustomFunction):
    """
    This class chains two Ramp functions
    """
    def __init__(self, first_p_ramp, second_p_ramp):
        """
        :param first_p_ramp: first ramp
        :type first_p_ramp: Ramp
        :param second_p_ramp: second ramp
        :type second_p_ramp: Ramp
        """
        self.__ramp_1 = first_p_ramp
        self.__ramp_2 = second_p_ramp

        if self.__ramp_1.end_time > self.__ramp_2.start_time:
            raise ValueError("""Cannot go into the past."""
                             """You have to build the ramp with increasing times""")

        if self.__ramp_1.second_value != self.__ramp_2.first_value:
            raise ValueError("""Please use a continuous law!""")

    def evaluate(self, time, *args, **kwargs):
        if time <= self.__ramp_2.start_time:
            return self.__ramp_1.evaluate(time)
        return self.__ramp_2.evaluate(time)
