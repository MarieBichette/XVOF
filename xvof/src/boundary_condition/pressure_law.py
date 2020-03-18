# -*- coding: utf-8 -*-
"""
Implements the PressureLaw abstract base class
"""
from abc import abstractmethod
from xvof.src.boundary_condition.boundary_condition import BoundaryCondition


class PressureLaw(BoundaryCondition):
    """
    An abstract base class for pressure law
    """
    def __init__(self):
        super(PressureLaw, self).__init__()
        self._type = "pressure"

    @abstractmethod
    def evaluate(self, time, *args, **kwargs):
        """
        Evaluate the pressure at the given time
        """
        pass
