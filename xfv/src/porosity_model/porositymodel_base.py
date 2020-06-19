# -*- coding: utf-8 -*-
"""
Interface for the shear modulus calculation
"""
from abc import abstractmethod
import numpy as np


class PorosityModelBase:  # pylint: disable=too-few-public-methods
    """
    Abstract class for the porosity model computation
    """
    def __init__(self, *args, **kwargs):
        """
        Build the porosity model
        """
        pass

    @abstractmethod
    def compute_porosity(self, delta_t: float, porosity: np.array,
                         pressure: np.array)->np.array:
        """
        Compute the new value of the porosity
        """
