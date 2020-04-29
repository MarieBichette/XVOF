# -*- coding: utf-8 -*-
"""
Interface for the shear modulus calculation
"""
from abc import abstractmethod
import numpy as np


class ShearModulus:  # pylint: disable=too-few-public-methods
    """
    Abstract class for the shear modulus computation
    """
    def __init__(self, initial_value):
        self.shear_modulus = initial_value

    @abstractmethod
    def compute(self, density: np.array) -> np.array:
        """
        Compute the new value of shear modulus
        """
