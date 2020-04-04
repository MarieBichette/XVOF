# -*- coding: iso-8859-1 -*-
"""
Definition of EquationOfStateBase interface
"""
from abc import ABCMeta, abstractmethod


class EquationOfStateBase(metaclass=ABCMeta):  # pylint: disable=too-few-public-methods
    """
    An interface for all equation of states
    """
    @abstractmethod
    def solve_volume_energy(self, specific_volume, internal_energy, pressure,  # pylint: disable=too-many-arguments
                            derivative, vson=None):
        """
        Solve the eos with [v, e] formulation
        """
