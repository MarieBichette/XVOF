"""
This module defines the classes that stores data read from the datafile and
needed to create UnloadingModel objects.
"""
from dataclasses import dataclass, field, asdict
from typing import Type

from xfv.src.cohesive_model_unloading.unloading_model_base import UnloadingModelBase
from xfv.src.cohesive_model_unloading.constant_stiffness_unloading import ConstantStiffnessUnloading
from xfv.src.cohesive_model_unloading.loss_of_stiffness_unloading import LossOfStiffnessUnloading


@dataclass
class UnloadingModelProps:  # pylint: disable=missing-class-docstring
    _unloading_model_class: Type[UnloadingModelBase] = field(init=False, repr=False)

    @staticmethod
    def dict_factory(obj):
        """
        Removes the classes (instance of type) that are inside obj
        """
        result = {}
        for key, value in obj:
            if not isinstance(value, type):
                result[key] = value
        return result

    def build_unloading_model_obj(self):
        """
        A factory that build and return the UnloadingModel object
        """
        return self._unloading_model_class(**asdict(self, dict_factory=self.dict_factory))


@dataclass  # pylint: disable=missing-class-docstring
class ConstantStiffnessUnloadingProps(UnloadingModelProps):
    slope: float
    _unloading_model_class = ConstantStiffnessUnloading


@dataclass  # pylint: disable=missing-class-docstring
class LossOfStiffnessUnloadingProps(UnloadingModelProps):
    _unloading_model_class = LossOfStiffnessUnloading
