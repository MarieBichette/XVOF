# -*- coding: utf-8 -*-
"""
Implements the SuccessivePressureRamp class
"""
from xvof.src.boundary_condition.pressure_law import PressureLaw


class SuccessivePressureRamp(PressureLaw):
    """
    This class defines a pressure boundary condition that chain two PressureRamp boundary conditions
    """
    def __init__(self, first_p_ramp, second_p_ramp):
        """
        :param first_p_ramp: first pressure ramp
        :type first_p_ramp: PressureRamp
        :param second_p_ramp: second pressure ramp
        :type second_p_ramp: PressureRamp
        """
        super(SuccessivePressureRamp, self).__init__()
        self.__ramp_1 = first_p_ramp
        self.__ramp_2 = second_p_ramp

        if self.__ramp_1.end_time > self.__ramp_2.start_time:
            raise ValueError("""Cannot go into the past."""
                             """You have to build the ramp with increasing times""")

        if self.__ramp_1.second_value != self.__ramp_2.first_value:
            raise ValueError("""Please use a continuous pressure law!""")

    def evaluate(self, time, *args, **kwargs):
        if time <= self.__ramp_2.start_time:
            return self.__ramp_1.evaluate(time)
        return self.__ramp_2.evaluate(time)
