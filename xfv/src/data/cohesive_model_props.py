"""
This module defines the classes that stores data read from the datafile and
needed to create CohesiveModel objects.
"""
from dataclasses import dataclass, field, asdict
from typing import Type

from xfv.src.data.unloading_model_props import UnloadingModelProps
from xfv.src.cohesive_model.cohesive_law_base import CohesiveZoneModelBase
from xfv.src.cohesive_model.linear_cohesive_law import LinearCohesiveZoneModel
from xfv.src.cohesive_model.bilinear_cohesive_law import BilinearCohesiveZoneModel
from xfv.src.cohesive_model.trilinear_cohesive_law import TrilinearCohesiveZoneModel


@dataclass  # pylint: disable=missing-class-docstring
class CohesiveModelProps:
    cohesive_strength: float
    critical_separation: float
    unloading_model: UnloadingModelProps
    _cohesive_model_class: Type[CohesiveZoneModelBase] = field(init=False, repr=False)

    @staticmethod
    def dict_factory(obj):
        """
        Removes the classes (instance of type) that are inside obj
        """
        result = {}
        for key, value in obj:
            if not isinstance(value, type):
                try:
                    value.build_unloading_model_obj()  # Case of unloading_model field
                except AttributeError:
                    pass
                result[key] = value
        return result

    def build_cohesive_model_obj(self):
        """
        Build and return the CohesiveModel object
        """
        return self._cohesive_model_class(**asdict(self, dict_factory=self.dict_factory))


@dataclass  # pylint: disable=missing-class-docstring
class LinearCohesiveZoneModelProps(CohesiveModelProps):
    _cohesive_model_class = LinearCohesiveZoneModel


@dataclass  # pylint: disable=missing-class-docstring
class BilinearCohesiveZoneModelProps(CohesiveModelProps):
    separation_1: float
    stress_1: float
    _cohesive_model_class = BilinearCohesiveZoneModel


@dataclass  # pylint: disable=missing-class-docstring
class TrilinearCohesiveZoneModelProps(CohesiveModelProps):
    separation_1: float
    stress_1: float
    separation_2: float
    stress_2: float
    _cohesive_model_class = TrilinearCohesiveZoneModel
