#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de base définissant une méthode d'incrémentation
"""
from abc import ABCMeta, abstractmethod


class NewtonRaphsonIncrementBase(object, metaclass=ABCMeta):
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
