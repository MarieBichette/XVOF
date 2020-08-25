# -*- coding: utf-8 -*-
"""
Implements the CustomFunction abstact base class
"""
from abc import ABCMeta, abstractmethod


class CustomFunction(metaclass=ABCMeta):
    """
    A generic CustomFunction class that should be used to derive
    more specific custom function class
    """

    __types = {"Pressure": [],
               "Velocity": []}

    @abstractmethod
    def evaluate(self, time, *args, **kwargs):  #pylint: disable=missing-docstring
        """
        Returns the value of the function evaluated at time

        :param time: the required time
        :param args: other arguments
        :param kwargs: other keywords arguments
        :return: the value
        """
        pass

    def register_pressure(self):
        """
        Register the current instance as a pressure
        """
        if self in CustomFunction.__types['Velocity']:
            raise ValueError("The instance {} is already registered as a velocity"
                             .format(self))
        if self not in CustomFunction.__types['Pressure']:
            CustomFunction.__types['Pressure'].append(self)

    def register_velocity(self):
        """
        Register the current instance as a velocity
        """
        if self in CustomFunction.__types['Pressure']:
            raise ValueError("The instance {} is already registered as a pressure"
                             .format(self))
        if self not in CustomFunction.__types['Velocity']:
            CustomFunction.__types['Velocity'].append(self)

    def is_pressure(self):
        """
        Return true if the instance is a pressure
        """
        return self in CustomFunction.__types['Pressure']

    def is_velocity(self):
        """
        Return true if the instance is a velocity
        """
        return self in CustomFunction.__types['Velocity']
