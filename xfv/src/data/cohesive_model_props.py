"""
This module defines the classes that stores data read from the datafile and
needed to create CohesiveModel objects.
"""
from dataclasses import dataclass, field, astuple
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
    __cohesive_model_class: Type[CohesiveZoneModelBase] = field(init=False, repr=False)

    @staticmethod
    def tuple_factory(obj):
        """
        Removes the classes (instance of type) that are inside obj
        """
        result = []
        for value in obj:
            if not isinstance(value, type):
                try:
                    value.build_unloading_model_obj()  # Case of unloading_model field
                except AttributeError:
                    pass
                result.append(value)
        return tuple(result)

    def build_cohesive_model_obj(self):
        """
        Build and return the CohesiveModel object
        """
        return self.build_cohesive_model_obj(*astuple(self, tuple_factory=self.tuple_factory))


@dataclass  # pylint: disable=missing-class-docstring
class LinearCohesiveZoneModelProps(CohesiveModelProps):
    __cohesive_model_class = LinearCohesiveZoneModel


@dataclass  # pylint: disable=missing-class-docstring
class BilinearCohesiveZoneModelProps(CohesiveModelProps):
    separation_1: float
    contrainte_1: float
    __cohesive_model_class = BilinearCohesiveZoneModel


@dataclass  # pylint: disable=missing-class-docstring
class TrilinearCohesiveZoneModelProps(CohesiveModelProps):
    separation_1: float
    contrainte_1: float
    separation_2: float
    contrainte_2: float
    __cohesive_model_class = TrilinearCohesiveZoneModel