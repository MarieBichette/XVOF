"""
This module defines the classes that stores data read from the datafile and
needed to create ShearModulus objects.
"""
from dataclasses import dataclass, field, asdict
from typing import Type

from xfv.src.rheology.shearmodulus import ShearModulus
from xfv.src.rheology.constantshearmodulus import ConstantShearModulus


@dataclass  # pylint: disable=missing-class-docstring
class ShearModulusProps:
    _shear_modulus_class: Type[ShearModulus] = field(init=False, repr=False)

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

    def build_shear_modulus_obj(self):
        """
        A factory that builds and returns the ShearModulus object
        """
        return self._shear_modulus_class(**asdict(self, dict_factory=self.dict_factory))


@dataclass  # pylint: disable=missing-class-docstring, too-many-instance-attributes
class ConstantShearModulusProps(ShearModulusProps):
    init_value: float
    _shear_modulus_class = ConstantShearModulus
