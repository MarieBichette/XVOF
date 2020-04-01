"""
This module defines the classes that stores data read from the datafile and
needed to create PlasticityCriterion objects.
"""
from dataclasses import dataclass, field
from typing import Type

from xfv.src.data.type_checked_dataclass import TypeCheckedDataClass
from xfv.src.plasticitycriterion.plasticitycriterion import PlasticityCriterion
from xfv.src.plasticitycriterion.vonmises import VonMisesCriterion


@dataclass  # pylint: disable=missing-class-docstring
class PlasticityCriterionProps(TypeCheckedDataClass):
    _plasticity_criterion_class: Type[PlasticityCriterion] = field(init=False, repr=False)

    def build_plasticity_criterion_obj(self):
        """
        A factory that builds the PlasticityCriterion object
        """
        return self._plasticity_criterion_class()


@dataclass  # pylint: disable=missing-class-docstring
class VonMisesCriterionProps(PlasticityCriterionProps):
    _plasticity_criterion_class = VonMisesCriterion
