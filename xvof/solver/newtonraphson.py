#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant un solveur non linéaire de type Newton Raphson

:todo: Mettre les critères de convergence dans le XML des données
"""
import numpy as np
from xvof.solver.incrementmethods.classicalnewtonraphson import ClassicalNewtonRaphsonIncrement
from xvof.solver.newtonraphsonbase import NewtonRaphsonBase

# EPSILON = 1.0e-09
# PRECISION = 1.0e-10
# EPSILON = 1.0e-08
# PRECISION = 1.0e-09
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
        # Variable du Newton
        var_i = init_variable
        var_iplus1 = np.zeros(var_i.shape, dtype=np.float64, order='C')
        # Critère de convergence
        not_conv = np.array([True for i in xrange(len(var_i))])
        # Nombre d'itérations
        nit = 0
        #
        while not_conv.any() and nit < self.nb_iterations_max:
            (func_i, dfunc_i_surde) = self.function.computeFunctionAndDerivative(var_i, not_conv)
            # Correction
            delta = self._increment_method.computeIncrement(func_i, dfunc_i_surde)
            var_iplus1[not_conv] = var_i[not_conv] + delta
            nit += 1
            not_conv[not_conv] = abs(func_i) >= EPSILON * abs(delta) + PRECISION
            if not not_conv.any():
                res = var_i
                break
            # Incrémentation
            var_i = var_iplus1
        if nit == self.nb_iterations_max:
            print "Erreur de convergence du NR"
            print "func_i=", func_i
            print "delta=", delta
            print "nit=", nit
            raise ValueError()
        return res
