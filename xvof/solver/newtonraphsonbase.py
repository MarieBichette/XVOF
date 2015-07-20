#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de base définissant un solveur non linéaire 1D de type Newton-Raphson
"""
from abc import ABCMeta, abstractmethod


class NewtonRaphsonBase(object):
    """
    Une classe de base pour les solveurs non linéaire 1D de type Newton-Raphson
    """
    # pylint: disable=abstract-class-not-used
    # Nécessaire pour spécifier l'interface
    __metaclass__ = ABCMeta
    #

    def __init__(self, function_to_vanish, variable, nb_iterations_max, increment_method):
        self.__function_to_vanish = function_to_vanish
        self.__variable_init = variable
        self.__nb_iterations_max = nb_iterations_max
        self._increment_method = increment_method

    @property
    def function(self):
        '''
        Renvoie la fonction à annuler
        '''
        return self.__function_to_vanish

    @property
    def init_variable(self):
        '''
        Retourne la variable initialisée
        '''
        return self.__variable_init

    @property
    def nb_iterations_max(self):
        '''
        Retourne le nombre maximal d'itérations autorisé
        '''
        return self.__nb_iterations_max

    @abstractmethod
    def computeSolution(self):
        '''
        Résolution du solveur
        '''
        raise NotImplementedError
