"""
This module defines the classes that stores data read from the datafile and
needed to create YieldStress objects.
"""
from dataclasses import dataclass, field, asdict
from typing import Type

from xfv.src.data.type_checked_dataclass import TypeCheckedDataClass
from xfv.src.rheology.yieldstress import YieldStress
from xfv.src.rheology.constantyieldstress import ConstantYieldStress


@dataclass  # pylint: disable=missing-class-docstring
class YieldStressProps(TypeCheckedDataClass):
    _yield_stress_class: Type[YieldStress] = field(init=False, repr=False)

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

    def build_yield_stress_obj(self):
        """
        A factory that builds and returns the YieldStress object
        """
        return self._yield_stress_class(**asdict(self, dict_factory=self.dict_factory))


@dataclass  # pylint: disable=missing-class-docstring, too-many-instance-attributes
class ConstantYieldStressProps(YieldStressProps):
    init_value: float
    _yield_stress_class = ConstantYieldStress
