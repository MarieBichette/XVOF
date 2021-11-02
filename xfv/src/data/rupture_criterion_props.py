"""
This module defines the classes that stores data read from the datafile and
needed to create RuptureCriterion objects.
"""
from dataclasses import dataclass, field, asdict
from typing import Type, Optional

from xfv.src.data.type_checked_dataclass import TypeCheckedDataClass
from xfv.src.rupturecriterion.rupturecriterion import RuptureCriterion
from xfv.src.rupturecriterion.damage_criterion import DamageCriterion
from xfv.src.rupturecriterion.porosity_criterion import PorosityCriterion
from xfv.src.rupturecriterion.halfrodcomparison import HalfRodComparisonCriterion
from xfv.src.rupturecriterion.maximalstress import MaximalStressCriterion
from xfv.src.rupturecriterion.minimumpressure import MinimumPressureCriterion
from xfv.src.rupturecriterion.nonlocalstress import NonLocalStressCriterion
from xfv.src.rupturecriterion.nonlocalstresswithmax import NonLocalStressCriterionWithMax


@dataclass  # pylint: disable=missing-class-docstring
class RuptureCriterionProps(TypeCheckedDataClass):
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
    sigma_max: Optional[float]  # optional to personalize the error message
    _rupture_criterion_class = MaximalStressCriterion

    def __post_init__(self):
        super().__post_init__()  #  typecheck first
        self._ensure_defined('sigma_max', 'MaximalStressCriterionProps',
                             'failure/failure-criterion/value')  # ensures that exists


@dataclass  # pylint: disable=missing-class-docstring
class DamageCriterionProps(RuptureCriterionProps):
    d_limite: Optional[float]  # optional to personalize the error message
    _rupture_criterion_class = DamageCriterion

    def __post_init__(self):
        super().__post_init__()  # typecheck first
        self._ensure_defined('d_limite', 'DamageCriterionProps',
                             'failure/failure-criterion/value')


@dataclass  # pylint: disable=missing-class-docstring
class PorosityCriterionProps(RuptureCriterionProps):
    p_limit: Optional[float]  # optional to personalize the error message
    _rupture_criterion_class = PorosityCriterion

    def __post_init__(self):
        super().__post_init__()  # typecheck first
        self._ensure_defined('p_limit', 'PorosityCriterionProps',
                             'failure/failure-criterion/value')


@dataclass  # pylint: disable=missing-class-docstring
class HalfRodComparisonCriterionProps(RuptureCriterionProps):
    ruptured_cell_index: Optional[int]  # optional to personalize the error message
    _rupture_criterion_class = HalfRodComparisonCriterion

    def __post_init__(self):
        super().__post_init__()  # typecheck first
        self._ensure_defined('ruptured_cell_index', 'HalfRodComparisonCriterionProps',
                             'failure/failure-criterion/index')  # ensures that exists


@dataclass  # pylint: disable=missing-class-docstring
class MinimumPressureCriterionProps(RuptureCriterionProps):
    pmin: Optional[float]  # optional to personalize the error message
    _rupture_criterion_class = MinimumPressureCriterion

    def __post_init__(self):
        super().__post_init__()  # typecheck first
        self._ensure_defined('pmin', 'MinimumPressureCriterionProps',
                             'failure/failure-criterion/value')  # ensures that exists


@dataclass  # pylint: disable=missing-class-docstring
class NonLocalStressCriterionProps(RuptureCriterionProps):
    value: Optional[float]  # optional to personalize the error message
    radius: Optional[float]
    _rupture_criterion_class = NonLocalStressCriterion

    def __post_init__(self):
        super().__post_init__()  # typecheck first
        self._ensure_defined('value', 'NonLocalStressCriterionProps',
                             'failure/failure-criterion/value')  # ensures that exists

@dataclass  # pylint: disable=missing-class-docstring
class NonLocalStressCriterionWithMaxProps(RuptureCriterionProps):
    value: Optional[float]  # optional to personalize the error message
    radius: Optional[float]
    _rupture_criterion_class = NonLocalStressCriterionWithMax

    def __post_init__(self):
        super().__post_init__()  # typecheck first
        self._ensure_defined('value', 'NonLocalStressCriterionWithMaxProps',
                             'failure/failure-criterion/value')  # ensures that exists