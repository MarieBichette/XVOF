"""
This module defines the classes that stores data read from the datafile and
needed to create UnloadingModel objects.
"""
from dataclasses import dataclass, field
from typing import Type

from xfv.src.cohesive_model.unloading_model_base import UnloadingModelBase
from xfv.src.cohesive_model.zero_force_unloading_model import ZeroForceUnloadingModel
from xfv.src.cohesive_model.progressive_unloading_model import ProgressiveUnloadingModel
from xfv.src.cohesive_model.loss_of_stiffness_unloading_model import LossOfStiffnessUnloadingModel


@dataclass
class UnloadingModelProps:  # pylint: disable=missing-class-docstring
    slope: float
    cohesive_strength: float
    __unloading_model_class: Type[UnloadingModelBase] = field(init=False, repr=False)

    def build_unloading_model_obj(self):
        """
        A factory that build and return the UnloadingModel object
        """
        return self.__unloading_model_class(self.slope, self.cohesive_strength)


@dataclass  # pylint: disable=missing-class-docstring
class ZeroForceUnloadingModelProps(UnloadingModelProps):
    __unloading_model_class = ZeroForceUnloadingModel


@dataclass  # pylint: disable=missing-class-docstring
class ProgressiveUnloadingModelProps(UnloadingModelProps):
    __unloading_model_class = ProgressiveUnloadingModel


@dataclass  # pylint: disable=missing-class-docstring
class LossOfStiffnessUnloadingModelProps(UnloadingModelProps):
    __unloading_model_class = LossOfStiffnessUnloadingModel
