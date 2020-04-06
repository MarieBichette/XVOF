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
        var_i = np.ndarray(init_variable.shape, dtype=np.float64, order='C')
        var_i[:] = init_variable
        func_i = np.ndarray(var_i.shape, dtype=np.float64, order='C')
        dfunc_i_surde = np.ndarray(var_i.shape, dtype=np.float64, order='C')

        # Newton's parameters initialization
        # Convergence criterion cell by cell
        above_crit = np.ndarray(var_i.shape, dtype=bool, order='C')
        above_crit[:] = True
        is_conv = False
        nit = 0  # Number of iterations

        while not is_conv and nit < self.nb_iterations_max:
            func_i[above_crit], dfunc_i_surde[above_crit] = (
                self.function.computeFunctionAndDerivative(var_i, above_crit))
            delta = self._increment_method.computeIncrement(func_i, dfunc_i_surde)
            above_crit = abs(func_i) >= EPSILON * abs(delta) + PRECISION
            is_conv = not above_crit.any()
            if is_conv:
                break
            var_i += delta
            nit += 1

        # Error if non convergence
        if nit == self.nb_iterations_max:
            msg = ("Erreur de convergence du NR\n"
                   f"func_i = {func_i}\n"
                   f"delta = {delta}\n"
                   f"nit = {nit}")
            raise ValueError(msg)
        return var_i
