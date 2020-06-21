# -*- coding: utf-8 -*-
"""
Class for the computation of yield stress
"""
from abc import abstractmethod
import numpy as np


class YieldStress:  # pylint: disable=too-few-public-methods
    """
    Interface for yield stress computation
    """
    def __init__(self, initial_value):
        """
        Initialization of the class
        :param initial_value: yield_stress initial value
        """
        self.yield_stress = initial_value

    @abstractmethod
    def compute(self, density: np.array) -> np.array:
        """
        Compute the new value of shear modulus

        :param density: the current density
        :return: the computed yield stress
        """
