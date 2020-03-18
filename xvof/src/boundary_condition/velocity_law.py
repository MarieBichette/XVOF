# -*- coding: utf-8 -*-
"""
Implements the VelocityLaw abstract base class
"""
from abc import abstractmethod
from xvof.src.boundary_condition.boundary_condition import BoundaryCondition


class VelocityLaw(BoundaryCondition):
    """
    An abstract base class for pressure law
    """
    def __init__(self):
        super(VelocityLaw, self).__init__()
        self._type = "velocity"

    @abstractmethod
    def evaluate(self, time, *args, **kwargs):
        """
        Evaluate the velocity at the given time
        """
        pass
