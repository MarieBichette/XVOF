"""
This module defines the classes that stores data read from the datafile and
needed to create CohesiveModel objects.
"""
from dataclasses import dataclass

import numpy as np

from xfv.src.data.type_checked_dataclass import TypeCheckedDataClass
from xfv.src.data.unloading_model_props import UnloadingModelProps
from xfv.src.cohesive_model.cohesive_zone_model import CohesiveZoneModel


@dataclass  # pylint: disable=missing-class-docstring
class CohesiveZoneModelProps(TypeCheckedDataClass):
    cohesive_strength: float
    critical_opening: float
    cohesive_zone_model_name: str
    unloading_model: UnloadingModelProps
    dissipated_energy: float
    purcentage_internal_energy: float
    _cohesive_zone_model_class = CohesiveZoneModel
    
    def _build_cohesive_law(self):
        """
        Build the cohesive law that is needed by the CohesiveModel
        """
        raise NotImplementedError("This is an abstract method!")
    
    def cohesive_model_name(self):
        return self.cohesive_zone_model_name

    def purcentage(self):
        return self.purcentage_internal_energy
    
    def dissip_energy(self):
        return self.dissipated_energy

    def build_cohesive_model_obj(self):
        """
        Build and return the CohesiveModel object
        """
        return self._cohesive_zone_model_class(self._build_cohesive_law(),
                                               self.unloading_model.build_unloading_model_obj(),
                                               self.cohesive_model_name(), self.dissip_energy(),
                                               self.purcentage())

@dataclass  # pylint: disable=missing-class-docstring
class LinearMixedCohesiveZoneModelProps(CohesiveZoneModelProps):

    def _build_cohesive_law(self):
        return np.array([
            [0, self.cohesive_strength],
            [self.critical_opening, 0]])


@dataclass  # pylint: disable=missing-class-docstring  
class LinearCohesiveZoneModelProps(CohesiveZoneModelProps):

    def _build_cohesive_law(self):
        return np.array([
            [0, self.cohesive_strength],
            [self.critical_opening, 0]])


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
