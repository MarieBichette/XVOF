# -*- coding: utf-8 -*-
# pylint: disable=too-few-public-methods
"""
Class définissant la correction amortie appliquée sur la variable d'un Newton Raphson
"""
from xfv.src.solver.incrementmethods.newtonraphsonincrementbase import NewtonRaphsonIncrementBase


class DampedNewtonRaphsonIncrement(NewtonRaphsonIncrementBase):
    """
    Class définissant un incrément amorti de l'algorithme de Newton-Raphson
    """
    def __init__(self, damping_coefficient=0.9):
        super().__init__()
        self._damping_coefficient = damping_coefficient

    def computeIncrement(self, function_value, derivative_function_value):
        """
        Increment computation
        """
        return - self._damping_coefficient * function_value / derivative_function_value
