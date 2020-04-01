"""
This module defines the classes that stores data read from the datafile and
needed to create RuptureCriterion objects.
"""
from dataclasses import dataclass, field, asdict
from typing import Type

from xfv.src.rupturecriterion.rupturecriterion import RuptureCriterion
from xfv.src.rupturecriterion.damage_criterion import DamageCriterion
from xfv.src.rupturecriterion.halfrodcomparison import HalfRodComparisonCriterion
from xfv.src.rupturecriterion.maximalstress import MaximalStressCriterion
from xfv.src.rupturecriterion.minimumpressure import MinimumPressureCriterion


@dataclass  # pylint: disable=missing-class-docstring
class RuptureCriterionProps:
    _rupture_criterion_class: Type[RuptureCriterion] = field(init=False, repr=False)

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

    def build_rupture_criterion_obj(self):
        """
        Build and return the CohesiveModel object
        """
        return self._rupture_criterion_class(**asdict(self, dict_factory=self.dict_factory))


@dataclass  # pylint: disable=missing-class-docstring
class MaximalStressCriterionProps(RuptureCriterionProps):
    sigma_max: float
    _rupture_criterion_class = MaximalStressCriterion


@dataclass  # pylint: disable=missing-class-docstring
class DamageCriterionProps(RuptureCriterionProps):
    d_limite: float
    _rupture_criterion_class = DamageCriterion


@dataclass  # pylint: disable=missing-class-docstring
class HalfRodComparisonCriterionProps(RuptureCriterionProps):
    ruptured_cell_index: int
    _rupture_criterion_class = HalfRodComparisonCriterion


@dataclass  # pylint: disable=missing-class-docstring
class MinimumPressureCriterionProps(RuptureCriterionProps):
    pmin: float
    _rupture_criterion_class = MinimumPressureCriterion
