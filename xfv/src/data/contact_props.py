"""
This module defines the classes that stores data read from the datafile and
needed to create Contact objects.
"""
from dataclasses import dataclass, field, asdict
from typing import Type

from xfv.src.data.type_checked_dataclass import TypeCheckedDataClass
from xfv.src.contact.contact_base import ContactBase
from xfv.src.contact.penalty import PenaltyContact
from xfv.src.contact.lagrange_multiplier import LagrangianMultiplierContact


@dataclass  # pylint: disable=missing-class-docstring
class ContactProps(TypeCheckedDataClass):
    _contact_class: Type[ContactBase] = field(init=False, repr=False)

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

    def build_contact_obj(self):
        """
        A factory that builds and returns the YieldStress object
        """
        return self._contact_class(**asdict(self, dict_factory=self.dict_factory))


@dataclass  # pylint: disable=missing-class-docstring, too-many-instance-attributes
class PenaltyContactProps(ContactProps):
    penalty_stiffness: float
    _contact_class = PenaltyContact


@dataclass  # pylint: disable=missing-class-docstring, too-many-instance-attributes
class LagrangianMultiplierProps(ContactProps):
    _contact_class = LagrangianMultiplierContact