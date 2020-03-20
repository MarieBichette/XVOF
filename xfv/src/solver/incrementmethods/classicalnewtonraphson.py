#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant la correction classique appliquée sur la variable d'un Newton Raphson
"""
from xvof.src.solver.incrementmethods.newtonraphsonincrementbase import NewtonRaphsonIncrementBase


class ClassicalNewtonRaphsonIncrement(NewtonRaphsonIncrementBase):
    """
    Classe définissant un incrément classique de l'algorithme de Newton-Raphson
    """
    def __init__(self):
        super(ClassicalNewtonRaphsonIncrement, self).__init__()

    def computeIncrement(self, function_value, derivative_function_value):
        """
        Increment computation
        """
        return - function_value / derivative_function_value
