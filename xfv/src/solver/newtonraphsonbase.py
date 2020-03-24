# -*- coding: utf-8 -*-
"""
Implements the NewtonRaphsonBase abstract base class
"""
from abc import ABCMeta, abstractmethod


class NewtonRaphsonBase(metaclass=ABCMeta):
    """
    This a base class for 1D non linear Newton-Raphson solver
    """
    def __init__(self, function_to_vanish, nb_iterations_max, increment_method):
        self.__function_to_vanish = function_to_vanish
        self.__nb_iterations_max = nb_iterations_max
        self._increment_method = increment_method

    @property
    def function(self):
        """
        Returns the function to vanish
        """
        return self.__function_to_vanish

    @property
    def nb_iterations_max(self):
        """
        Returns the maximum iteration number
        """
        return self.__nb_iterations_max

    @abstractmethod
    def compute_solution(self, init_variable):
        """
        Compute the solution
        """
        raise NotImplementedError
