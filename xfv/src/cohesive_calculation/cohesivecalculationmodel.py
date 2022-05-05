# -*- coding: utf-8 -*-
"""
Interface for the cohesive calculation model
"""
from abc import abstractmethod
import numpy as np


class CohesiveCalculationModel:  # pylint: disable=too-few-public-methods
    """
    Abstract class for the cohesive calculation model
    """


    @abstractmethod
    def get_values(energy, stress, mass, section, ind, cohesive_model)-> float:
        """
        Compute and return critical strength and critical separation

        :param energy: array for internal energy of cells [J/kg]
        :param stress: array for stress of cells [Pa]
        :param mass: array for mass of cells [kg]
        :param section: float for section of material [m^2]
        :param ind: integer for indice of cell
        :cohesive_model: type of cohesive model 
        """
