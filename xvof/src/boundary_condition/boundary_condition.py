# -*- coding: utf-8 -*-
"""
Implements the BoundaryCondition abstact base class
"""
from abc import ABCMeta, abstractmethod


class BoundaryCondition(object):
    """
    A generic BoundaryCondition class that should be used to derive
    more specific boundary condition class
    """
    __metaclass__ = ABCMeta

    __types = {"Pressure": [],
               "Velocity": []}

    @abstractmethod
    def evaluate(self, time, *args, **kwargs):  #pylint: disable=missing-docstring
        pass

    @classmethod
    def register_pressure_bc(cls, bc_obj):
        """
        Register the bc_obj as a pressure boundary condition
        """
        cls.__types['Pressure'].append(bc_obj)

    @classmethod
    def register_velocity_bc(cls, bc_obj):
        """
        Register the bc_obj as a velocity boundary condition
        """
        cls.__types['Velocity'].append(bc_obj)

    def is_pressure_boundary_condition(self):
        """
        Return true if the instance is a pressure type boundary conditions
        """
        return self in BoundaryCondition.__types['Pressure']

    def is_velocity_boundary_condition(self):
        """
        Return true if the instance is a velocity type boundary conditions
        """
        return self in BoundaryCondition.__types['Velocity']
