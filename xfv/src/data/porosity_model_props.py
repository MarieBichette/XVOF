"""
This module defines the classes that stores data read from the datafile and
needed to create JohnsonModel objects.
"""
from abc import abstractmethod
from dataclasses import dataclass, field, asdict
from typing import Type
import numpy as np

from xfv.src.data.type_checked_dataclass import TypeCheckedDataClass
from xfv.src.porosity_model.porositymodel_base import PorosityModelBase
from xfv.src.porosity_model.johnsonmodel import JohnsonModel

@dataclass  # pylint: disable=missing-class-docstring
class PorosityModelProps(TypeCheckedDataClass):
    _porosity_model_class: Type[PorosityModelBase] = field(init=False, repr=False)

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

    def build_porosity_model_obj(self):
        """
        A factory that builds and returns the PorosityModel object
        """
        return self._porosity_model_class(**asdict(self, dict_factory=self.dict_factory))

    @abstractmethod
    def compute_porosity(self, delta_t: float, porosity: np.array,
                                               pressure: np.array) -> np.array:
        """
        Compute the new value of porosity
        """

@dataclass  # pylint: disable=missing-class-docstring
class JohnsonModelProps(PorosityModelProps):
    initial_porosity_for_johnson: float
    effective_strength_for_johnson: float
    viscosity_for_johnson: float  # pylint: disable=invalid-name
    _porosity_model_class = JohnsonModel

    def __post_init__(self):
        super().__post_init__()  # typecheck first
        self._ensure_strict_positivity('effective_strength_for_johnson', 'viscosity_for_johnson')

