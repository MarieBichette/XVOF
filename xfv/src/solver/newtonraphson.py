#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant un solveur non linéaire de type Newton Raphson

:todo: Mettre les critères de convergence dans le XML des données
"""
import numpy as np
from xfv.src.solver.incrementmethods.classicalnewtonraphson import ClassicalNewtonRaphsonIncrement
from xfv.src.solver.newtonraphsonbase import NewtonRaphsonBase

EPSILON = 1.0e-06
PRECISION = 1.0e-08  # valeurs par défaut dans A1 et dans A3

class NewtonRaphson(NewtonRaphsonBase):
    """
    Solveur non linéaire de type Newton Raphson
    """
    def __init__(self, function_to_vanish):
        super(NewtonRaphson, self).__init__(function_to_vanish, 100, ClassicalNewtonRaphsonIncrement())

    def setIncrementMethod(self, increment_method_obj):
        """
        Permet un changement de méthode d'incrémentation
        """
        self._increment_method = increment_method_obj

    def computeSolution(self, init_variable):
        """
        Algorithme de Newton-Raphson
        """
        if init_variable.size == 0:  # this case should never append but has been discovered in Unittests...
            msg = "Initialization variable for Newton has null size. Impossible to start Newton procedure."
            raise ValueError(msg)

        # Variable du Newton
        var_i = init_variable
        var_iplus1 = np.zeros(var_i.shape, dtype=np.float64, order='C')

        # Newton's parameters initialization
        not_conv = np.array([True for i in xrange(len(var_i))])  # Convergence criterion cell by cell
        nit = 0  # Number of iterations
        res = 0.  # result
        delta = 0.  # increment
        func_i = 0.  # function evaluation

        while not_conv.any() and nit < self.nb_iterations_max:
            (func_i, dfunc_i_surde) = self.function.computeFunctionAndDerivative(var_i, not_conv)
            # Correction
            delta = self._increment_method.computeIncrement(func_i, dfunc_i_surde)
            var_iplus1[not_conv] = var_i[not_conv] + delta
            not_conv[not_conv] = abs(func_i) >= EPSILON * abs(delta) + PRECISION
            if not not_conv.any():
                res = var_i
                break
            # Increment
            var_i = var_iplus1
            nit += 1

        # Error if non convergence
        if nit == self.nb_iterations_max:
            msg = "Erreur de convergence du NR"
            msg += "func_i=", func_i
            msg += "delta=", delta
            msg += "nit=", nit
            raise ValueError(msg)
        return res
