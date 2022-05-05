# -*- coding: utf-8 -*-
"""
Implementation of the LinearPercent class
"""

from distutils.errors import DistutilsInternalError
import numpy as np
from xfv.src.cohesive_calculation.cohesivecalculationmodel import CohesiveCalculationModel

class LinearPercent(CohesiveCalculationModel) :  # pylint: disable=too-few-public-methods
    """
    Class for cohesive linear Percent 
    """

    def get_values(energy, stress, mass, section, ind, cohesive_model) -> float :
        """
        Compute and return critical strength and critical separation

        :param energy: array for internal energy of cells [J/kg]
        :param stress: array for stress of cells [Pa]
        :param mass: array for mass of cells [kg]
        :param section: float for section of material [m^2]
        :param ind: integer for indice of cell
        :cohesive_model: type of cohesive model 
        """

        critical_strength = abs(stress[ind, 0])
        dissipated_energy = energy[ind]*mass[ind]*cohesive_model.purcentage
        critical_separation = 2.*dissipated_energy/(critical_strength*section)
        print('critical strength =', critical_strength)
        print('critical separation =', critical_separation)
        print('cohesive dissipated energy =', dissipated_energy)
        return critical_strength, critical_separation
