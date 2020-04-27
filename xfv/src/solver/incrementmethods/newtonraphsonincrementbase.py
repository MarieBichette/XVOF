# -*- coding: utf-8 -*-
"""
Classe de base définissant une méthode d'incrémentation
"""
from abc import ABCMeta, abstractmethod


class NewtonRaphsonIncrementBase(metaclass=ABCMeta):
    """
    Classe base définissant une méthode d'incrémentation
    """

    def __init__(self):
        pass

    @abstractmethod
    def computeIncrement(self, function_value, derivative_function_value):
        """
        Calcul de l'incrément
        """
        raise NotImplementedError
