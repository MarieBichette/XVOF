#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Interface définissant le contrat à remplir pur qu'ne fonction sont résolue par le solveur
non linéaire de tpye Newton Raphson
"""
from abc import ABCMeta, abstractmethod


class FunctionToSolveBase(object):
    """
    Interface définissant le contrat à remplir pur qu'ne fonction sont résolue par le solveur
    non linéaire de tpye Newton Raphson
    """
    # pylint: disable=abstract-class-not-used
    # Nécessaire pour spécifier l'interface
    __metaclass__ = ABCMeta

    def __init__(self):
        self._variables = None

    def setVariables(self, variables):
        '''
        Fixation de la valeur des variables
        '''
        if self._variables is None:
            self._variables = variables
        else:
            raise ValueError("Impossible de fixer deux fois les variables de la fonction à annuler!")

    def eraseVariables(self):
        '''
        Remise à None de _variables
        '''
        self._variables = None

    @abstractmethod
    def computeFunctionAndDerivative(self, variable_value):
        '''
        Renvoie la valeur de la fonction et sa dérivée
        '''
        raise NotImplementedError
