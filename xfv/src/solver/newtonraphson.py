# -*- coding: utf-8 -*-
"""
Implements the NewtonRaphson class

:todo: Put criterion convergence in the XML data
"""
import numpy as np
from xfv.src.solver.incrementmethods.classicalnewtonraphson import ClassicalNewtonRaphsonIncrement
from xfv.src.solver.newtonraphsonbase import NewtonRaphsonBase

EPSILON = 1.0e-06
PRECISION = 1.0e-08

class NewtonRaphson(NewtonRaphsonBase):
    """
    This class implements a Newton Raphson type non linear solver
    """
    def __init__(self, function_to_vanish):
        super(NewtonRaphson, self).__init__(
            function_to_vanish, 100, ClassicalNewtonRaphsonIncrement())

    def set_increment_method(self, increment_method_obj):
        """
        Allow a change of incrementation method
        """
        self._increment_method = increment_method_obj

    def compute_solution(self, init_variable):
        """
        Compute the solution through Newton-Raphson algorithm
        """
        # This case should never append but has been discovered in Unittests...
        if init_variable.size == 0:
            msg = ("Initialization variable for Newton has null size. "
                   "Impossible to start Newton procedure.")
            raise ValueError(msg)

        # Newton's variable
        var_i = init_variable
        var_iplus1 = np.ndarray(var_i.shape, dtype=np.float64, order='C')

        # Newton's parameters initialization
        # Convergence criterion cell by cell
        not_conv = np.ndarray(var_i.shape, dtype=bool, order='C')
        not_conv[:] = True
        is_conv = False
        nit = 0  # Number of iterations
        res = 0.  # result
        delta = 0.  # increment
        func_i = 0.  # function evaluation

        while not is_conv and nit < self.nb_iterations_max:
            func_i, dfunc_i_surde = self.function.computeFunctionAndDerivative(var_i, not_conv)
            # Correction
            delta = self._increment_method.computeIncrement(func_i, dfunc_i_surde)
            var_iplus1[not_conv] = var_i[not_conv] + delta
            not_conv[not_conv] = abs(func_i) >= EPSILON * abs(delta) + PRECISION
            is_conv = not not_conv.any()
            if is_conv:
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
