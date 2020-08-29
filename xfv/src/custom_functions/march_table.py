# -*- coding: utf-8 -*-
"""
Implements the MarchTable class
"""
import numpy as np
from xfv.src.custom_functions.custom_function import CustomFunction


class MarchTable(CustomFunction):
    """
    This class defines a value interpolated from a text file
    """
    def __init__(self, data_file):
        self.__data_file = data_file
        self.__data = np.loadtxt(self.__data_file)
        self.time_data = self.__data[:, 0]
        self.pressure_data = self.__data[:, 1]

    def evaluate(self, time, *args, **kwargs):
        """
        Return the value of the for the time given in argument

        :param time: current time
        :return: the value
        """
        if time > self.time_data[-1]:
            pressure = 0
        else:
            indice_limite_inf = np.where(self.time_data <= time)[0][-1]
            indice_limite_sup = np.where(self.time_data >= time)[0][0]
            if indice_limite_inf == indice_limite_sup:
                # current time is inside the table
                pressure = self.pressure_data[indice_limite_inf]
            else:
                # linear interpolation between the two nearest times present in the table
                delta_pressure = (self.pressure_data[indice_limite_sup] -
                                  self.pressure_data[indice_limite_inf])
                delta_time = self.time_data[indice_limite_sup] - self.time_data[indice_limite_inf]
                pente = delta_pressure / delta_time
                pressure = (self.pressure_data[indice_limite_inf] +
                            (time - self.time_data[indice_limite_inf]) * pente)
        return pressure
