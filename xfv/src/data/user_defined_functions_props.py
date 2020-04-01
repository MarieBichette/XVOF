"""
This module defines the classes that stores data read from the datafile and
needed to create CustomFunction objects.
"""
from abc import ABCMeta
from dataclasses import dataclass, field, astuple
from typing import Union, Type

from xfv.src.data.type_checked_dataclass import TypeCheckedDataClass
from xfv.src.custom_functions.custom_function import CustomFunction
from xfv.src.custom_functions.constant_value import ConstantValue
from xfv.src.custom_functions.ramp import Ramp
from xfv.src.custom_functions.two_steps import TwoSteps
from xfv.src.custom_functions.successive_ramp import SuccessiveRamp
from xfv.src.custom_functions.march_table import MarchTable


@dataclass
class UserDefinedFunctionProps(TypeCheckedDataClass, metaclass=ABCMeta):
    """
    This class defines the base class of all user defined function datas.
    A user defined function data class stores the data read from the datafile
    and needed to build the corresponding CustomFunction object
    """
    # Two choices here:
    # - Dot not annotate the type here because this variable should not be a field so that
    #   it doesn't appears through the use of astuple function
    #   => introduces confusion for mypy
    # - Or use type annotation in order to get mypy happy and customize the use of astuple
    #   function by using a function (tuple_factory) that remove all classes (instance of type)
    #   that are present in the fields of the dataclass
    #    => solution adopted here
    _custom_func_class: Type[CustomFunction] = field(init=False, repr=False)

    @staticmethod
    def tuple_factory(obj):
        """
        Removes the classes (instance of type) that are inside obj
        """
        result = []
        for value in obj:
            if not isinstance(value, type):
                result.append(value)
        return tuple(result)

    def build_custom_func(self) -> CustomFunction:
        """
        A factory that returns the CustomFunction object corresponding
        to the user defined function data
        """
        return self._custom_func_class(*astuple(self, tuple_factory=self.tuple_factory))


@dataclass  # pylint: disable=missing-class-docstring
class ConstantValueFunctionProps(UserDefinedFunctionProps):
    value: float
    _custom_func_class = ConstantValue


@dataclass  # pylint: disable=missing-class-docstring
class RampFunctionProps(UserDefinedFunctionProps):
    start_value: float
    end_value: float
    start_time: float
    end_time: float
    _custom_func_class = Ramp


@dataclass  # pylint: disable=missing-class-docstring
class TwoStepsFunctionProps(UserDefinedFunctionProps):
    first_value: float
    second_value: float
    critical_time: float
    _custom_func_class = TwoSteps


@dataclass  # pylint: disable=missing-class-docstring
class MarchTableFunctionProps(UserDefinedFunctionProps):
    file: str
    _custom_func_class = MarchTable


@dataclass  # pylint: disable=missing-class-docstring
class SuccessiveRampFunctionProps(UserDefinedFunctionProps):
    first_ramp: RampFunctionProps
    second_ramp: RampFunctionProps

    def build_custom_func(self) -> CustomFunction:
        first = self.first_ramp.build_custom_func()
        second = self.second_ramp.build_custom_func()
        return SuccessiveRamp(first, second)


UserDefinedFunctionPropsType = Union[ConstantValueFunctionProps, MarchTableFunctionProps,
                                     RampFunctionProps, SuccessiveRampFunctionProps,
                                     TwoStepsFunctionProps]
