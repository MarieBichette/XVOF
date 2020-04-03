"""
This module defines the classes that stores data read from the datafile and
needed to create EquationOfState objects.
"""
from dataclasses import dataclass, field, asdict
from typing import Type, Optional

from xfv.src.data.type_checked_dataclass import TypeCheckedDataClass
from xfv.src.equationsofstate.equationofstatebase import EquationOfStateBase
from xfv.src.equationsofstate.miegruneisen import MieGruneisen


@dataclass  # pylint: disable=missing-class-docstring
class EquationOfStateProps(TypeCheckedDataClass):
    _eos_class: Type[EquationOfStateBase] = field(init=False, repr=False)
    _eos_inst: Optional[EquationOfStateBase] = field(init=False, repr=False)

    @staticmethod
    def dict_factory(obj):
        """
        Removes the classes (instance of type) that are inside obj
        """
        result = {}
        for key, value in obj:
            if not key in ('_eos_class', '_eos_inst'):
                result[key] = value
        return result

    def build_eos_obj(self):
        """
        Build and return the CohesiveModel object
        """
        if self._eos_inst is None:
            self._eos_inst = self._eos_class(**asdict(self, dict_factory=self.dict_factory))
        return self._eos_inst

@dataclass  # pylint: disable=missing-class-docstring, too-many-instance-attributes
class MieGruneisenProps(EquationOfStateProps):
    czero: float
    S1: float  # pylint: disable=invalid-name
    S2: float  # pylint: disable=invalid-name
    S3: float  # pylint: disable=invalid-name
    rhozero: float
    grunzero: float
    b: float  # pylint: disable=invalid-name
    ezero: float
    _eos_class = MieGruneisen
    _eos_inst = None

    def __post_init__(self):
        super().__post_init__()  # typecheck first
        self._ensure_strict_positivity('czero', 'rhozero')
