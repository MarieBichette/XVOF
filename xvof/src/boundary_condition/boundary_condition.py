# -*- coding: utf-8 -*-
"""
Implements the BoundaryCondition abstact base class
"""
from abc import abstractmethod


class BoundaryCondition(object):
    """
    A generic BoundaryCondition class that should be used to derive
    more specific boundary condition class
    """
    def __init__(self):
        self._type = None

    def type(self):
        """
        Accessor on the boundary condition type (pressure or velocity)
        :return:
        """
        return self._type

    @abstractmethod
    def evaluate(self, time, *args, **kwargs):  #pylint: disable=missing-docstring
        pass
