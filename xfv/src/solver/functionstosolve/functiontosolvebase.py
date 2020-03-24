# -*- coding: iso-8859-1 -*-
"""
Implementing FunctionToSolveBase interface
"""
from abc import ABCMeta, abstractmethod


class FunctionToSolveBase(object, metaclass=ABCMeta):
    """
    An interface for functions that should be solve by Nexton-Raphson non linear solver
    """

    def __init__(self):
        self._variables = None

    def setVariables(self, variables):
        """
        Set of the value of the variables
        """
        if self._variables is None:
            self._variables = variables
        else:
            raise ValueError("Impossible de fixer deux fois les variables de la fonction à annuler!")

    def eraseVariables(self):
        """
        Reset of _variables
        """
        self._variables = None

    @abstractmethod
    def computeFunctionAndDerivative(self, var_value, mask):
        """
        Return the values of the function and its derivative

        :param var_value: value of the variable
        :param mask: boolean mask on which the function has to be evaluated
        """
        raise NotImplementedError
