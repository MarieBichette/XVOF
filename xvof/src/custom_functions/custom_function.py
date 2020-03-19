# -*- coding: utf-8 -*-
"""
Implements the CustomFunction abstact base class
"""
from abc import ABCMeta, abstractmethod


class CustomFunction(object):
    """
    A generic CustomFunction class that should be used to derive
    more specific custom function class
    """
    __metaclass__ = ABCMeta

    __types = {"Pressure": [],
               "Velocity": []}

    @abstractmethod
    def evaluate(self, time, *args, **kwargs):  #pylint: disable=missing-docstring
        pass

    def register_pressure(self):
        """
        Register the bc_obj as a pressure boundary condition
        """
        CustomFunction.__types['Pressure'].append(self)

    def register_velocity(self):
        """
        Register the bc_obj as a velocity boundary condition
        """
        CustomFunction.__types['Velocity'].append(self)

    def is_pressure_boundary_condition(self):
        """
        Return true if the instance is a pressure type boundary conditions
        """
        return self in CustomFunction.__types['Pressure']

    def is_velocity_boundary_condition(self):
        """
        Return true if the instance is a velocity type boundary conditions
        """
        return self in CustomFunction.__types['Velocity']
