"""
This module defines the classes that stores data read from the datafile and
needed to create CohesiveModel objects.
"""
from dataclasses import dataclass

import numpy as np

from xfv.src.data.type_checked_dataclass import TypeCheckedDataClass
from xfv.src.data.unloading_model_props import UnloadingModelProps
from xfv.src.cohesive_model.cohesive_zone_model import CohesiveZoneModel
from xfv.src.cohesive_calculation.lineardata import LinearData
from xfv.src.cohesive_calculation.linearenergy import LinearEnergy
from xfv.src.cohesive_calculation.linearpercent import LinearPercent
from xfv.src.cohesive_calculation.cohesivecalculationmodel import CohesiveCalculationModel

@dataclass  # pylint: disable=missing-class-docstring
class CohesiveZoneModelProps(TypeCheckedDataClass):
    cohesive_strength: float
    critical_opening: float
    _cohesive_zone_model_name: str
    unloading_model: UnloadingModelProps
    _dissipated_energy: float
    _purcentage: float
    _cohesive_zone_model_class = CohesiveZoneModel

    def _build_cohesive_law(self):
        """
        Build the cohesive law that is needed by the CohesiveModel
        """
        raise NotImplementedError("This is an abstract method!")


@dataclass  # pylint: disable=missing-class-docstring  
class LinearDataCohesiveZoneModelProps(CohesiveZoneModelProps):

    _cohesive_calculation_model = LinearData

    def _build_cohesive_law(self):
        """
        Build and return the CohesiveLaw
        """
        return np.array([
            [0, self.cohesive_strength],
            [self.critical_opening, 0]])   

    def build_cohesive_model_obj(self):
        """
        Build and return the CohesiveModel object
        """
        return self._cohesive_zone_model_class(self._build_cohesive_law(),
                                               self.unloading_model.build_unloading_model_obj(),
                                               self._cohesive_zone_model_name, self._dissipated_energy,
                                               self._purcentage, self._cohesive_calculation_model)

@dataclass  # pylint: disable=missing-class-docstring
class LinearPercentCohesiveZoneModelProps(CohesiveZoneModelProps):

    _cohesive_calculation_model = LinearPercent

    def _build_cohesive_law(self):
        return np.array([
            [0, self.cohesive_strength],
            [self.critical_opening, 0]])

    def build_cohesive_model_obj(self):
        """
        Build and return the CohesiveModel object
        """
        return self._cohesive_zone_model_class(self._build_cohesive_law(),
                                               self.unloading_model.build_unloading_model_obj(),
                                               self._cohesive_zone_model_name, self._dissipated_energy,
                                               self._purcentage, self._cohesive_calculation_model)


@dataclass  # pylint: disable=missing-class-docstring
class LinearEnergyCohesiveZoneModelProps(CohesiveZoneModelProps):

    _cohesive_calculation_model = LinearEnergy

    def _build_cohesive_law(self):
        return np.array([
            [0, self.cohesive_strength],
            [self.critical_opening, 0]])
    
    def build_cohesive_model_obj(self):
        """
        Build and return the CohesiveModel object
        """
        return self._cohesive_zone_model_class(self._build_cohesive_law(),
                                               self.unloading_model.build_unloading_model_obj(),
                                               self._cohesive_zone_model_name, self._dissipated_energy,
                                               self._purcentage,self._cohesive_calculation_model)


@dataclass  # pylint: disable=missing-class-docstring
class BilinearCohesiveZoneModelProps(CohesiveZoneModelProps):
    opening_1: float
    stress_1: float

    def _build_cohesive_law(self):
        return np.array([
            [0, self.cohesive_strength],
            [self.opening_1, self.stress_1],
            [self.critical_opening, 0]])


@dataclass  # pylint: disable=missing-class-docstring
class TrilinearCohesiveZoneModelProps(CohesiveZoneModelProps):
    opening_1: float
    stress_1: float
    opening_2: float
    stress_2: float
  
    def _build_cohesive_law(self):
        return np.array([
            [0, self.cohesive_strength],
            [self.opening_1, self.stress_1],
            [self.opening_2, self.stress_2],
            [self.critical_opening, 0]])
