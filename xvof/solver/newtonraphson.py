#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant un solveur non linéaire de type Newton Raphson
"""
from xvof.solver.newtonraphsonbase import NewtonRaphsonBase
from xvof.solver.incrementmethods.classicalnewtonraphson import ClassicalNewtonRaphsonIncrement


class NewtonRaphson(NewtonRaphsonBase):
    '''
    Solveur non linéaire de type Newton Raphson
    '''
    def __init__(self, function_to_vanish, variable):
        super(NewtonRaphson, self).__init__(function_to_vanish, variable, 10, ClassicalNewtonRaphsonIncrement())

    def setIncrementMethod(self, increment_method_obj):
        '''
        Permet un changement de méthode d'incrémentation
        '''
        self._increment_method = increment_method_obj

    def computeSolution(self):
        """
        Algorithme de Newton-Raphson
        """
        # Variable du Newton
        var_i = self.init_variable
        # Critère de convergence
        convergence = False
        # Nombre d'itérations
        nit = 0
        #
        while not convergence and (nit < self.nb_iterations_max):
            (func_i, dfunc_i_surde) = self.function.computeFunctionAndDerivative(var_i)
            # Correction
            var_iplus1 = var_i + self._increment_method.computeIncrement(func_i, dfunc_i_surde)
            nit += 1
            if abs(func_i) < 1e-09:
                convergence = True
                res = var_i
                break
            if abs(dfunc_i_surde) < 1.e-09:
                print "Sortie du NR par manque de pente :-)"
                convergence = True
                res = var_i
                break
            # Incrémentation
            var_i = var_iplus1
        if nit == self.nb_iterations_max:
            print "Erreur de convergence du NR"
            print "func_i=", func_i
            print "nit=", nit
            raise ValueError()
        return res
