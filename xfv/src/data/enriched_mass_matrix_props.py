"""
This module defines the classes that stores data read from the datafile and
needed to create EnrichedMassMatrix objects.
"""
from dataclasses import dataclass, field, asdict
from typing import Type

from xfv.src.data.type_checked_dataclass import TypeCheckedDataClass
from xfv.src.mass_matrix.enriched_mass_matrix import EnrichedMassMatrix
from xfv.src.mass_matrix.enriched_mass_matrix_consistent import EnrichedMassMatrixConsistent
from xfv.src.mass_matrix.enriched_mass_matrix_lump_menouillard \
    import EnrichedMassMatrixLumpMenouillard
from xfv.src.mass_matrix.enriched_mass_matrix_lump_sum import EnrichedMassMatrixLumpSum


@dataclass  # pylint: disable=missing-class-docstring
class EnrichedMassMatrixProps(TypeCheckedDataClass):
    _enr_mass_matrix_class: Type[EnrichedMassMatrix] = field(init=False, repr=False)

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

    def build_enriched_mass_matrix_obj(self):
        """
        A factory that builds and returns the YieldStress object
        """
        return self._enr_mass_matrix_class(**asdict(self, dict_factory=self.dict_factory))


@dataclass  # pylint: disable=missing-class-docstring, too-many-instance-attributes
class ConsistentMassMatrixProps(EnrichedMassMatrixProps):
    _enr_mass_matrix_class = EnrichedMassMatrixConsistent


@dataclass  # pylint: disable=missing-class-docstring, too-many-instance-attributes
class LumpMenouillardMassMatrixProps(EnrichedMassMatrixProps):
    _enr_mass_matrix_class = EnrichedMassMatrixLumpMenouillard


@dataclass  # pylint: disable=missing-class-docstring, too-many-instance-attributes
class LumpSumMassMatrixProps(EnrichedMassMatrixProps):
    _enr_mass_matrix_class = EnrichedMassMatrixLumpSum
