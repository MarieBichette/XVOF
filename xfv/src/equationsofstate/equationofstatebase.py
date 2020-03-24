# -*- coding: iso-8859-1 -*-
"""
Definition of EquationOfStateBase interface
"""
from abc import ABCMeta, abstractmethod


class EquationOfStateBase(object, metaclass=ABCMeta):
    """
    An interface for all equation of states
    """

    @abstractmethod
    def solveVolumeEnergy(self, specific_volume, internal_energy):
        """
        Solve the eos with [v, e] formulation
        """
        pass

    @abstractmethod
    def solveVolumeTemperature(self, specific_volume, temperature):
        """
        Solve the eos with [v, T] formulation
        """
        pass

    @abstractmethod
    def solveVolumePressure(self, specific_volume, pressure):
        """
        Solve the eos with [v, P] formulation
        """
        pass
