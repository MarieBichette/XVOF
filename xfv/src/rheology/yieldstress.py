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
    def __init__(self, init_value, init_shear_modulus, Y_max, beta, m):
        """
        Initialization of the constant yield stress class

        :param init_value: initial yield stress
        """
        self.init_value = init_value
        self.init_shear_modulus = init_shear_modulus
        self.Y_max = Y_max
        self.beta = beta
        self.m = m 

    @abstractmethod
    def compute(self, density: np.array, strain_plastic_eq, G) -> np.array:
        """
        Compute the new value of shear modulus

        :param density: the current density
        :return: the computed yield stress
        """
