#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Class définissant la correction amortie appliquée sur la variable d'un Newton Raphson
"""
from xvof.solver.incrementmethods.newtonraphsonincrementbase import NewtonRaphsonIncrementBase


class DampedNewtonRaphsonIncrement(NewtonRaphsonIncrementBase):
    """
    Class définissant un incrément amorti de l'algorithme de Newton-Raphson
    """
    def __init__(self, damping_coefficient=0.9):
        super(DampedNewtonRaphsonIncrement, self).__init__()
        self._damping_coefficient = damping_coefficient

    def computeIncrement(self, function_value, derivative_function_value):
        """
        Increment computation
        """
        return - self._damping_coefficient * function_value / derivative_function_value
