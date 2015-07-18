#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de base abstraite (interface) définissant une équation d'état
"""
from abc import ABCMeta, abstractmethod


class EquationOfStateBase(object):
    """
    Une interface pour les équations d'état
    """
    # Nécessaire pour spécifier l'interface
    __metaclass__ = ABCMeta
    #

    def __init__(self):
        pass

    @abstractmethod
    def solveVolumeEnergy(self, specific_volume, internal_energy):
        """
        Résolution en formulation V-E
        """
        pass

    @abstractmethod
    def solveVolumeTemperature(self, specific_volume, temperature):
        """
        Résolution en formulation V-T
        """
        pass

    @abstractmethod
    def solveVolumePressure(self, specific_volume, pressure):
        """
        Résolution en formulation V-P
        """
        pass
